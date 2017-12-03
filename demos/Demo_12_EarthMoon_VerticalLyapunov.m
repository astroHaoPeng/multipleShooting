%Demonstration of multiple-shooting method to generate vertical Lyapunov orbits in the Earth-Moon system.
%
% 1. generate 3rd order approximation of vertical Lyapunov orbit
% 2. shoot in CRTBP
% 3. refine in ephemeris model using multiple-shooting
%
%Could be used for Sun-Earth/Moon system by
% must use corresponding 'lengthUnit', 'timeUnit'
% must use corresponding frame kernels 'SE_EarthCenteredRotation.fk', 'EarthCenteredInertial.fk'


%%
addpath('../ephemeris','-begin');
addpath('../lyapunov','-begin');
addpath('../rtbp','-begin');
addpath('../common','-begin');
addpath('../','-begin');
clear;
% parpool;

%% define constants

% G * mass
cspice_kclear;
cspice_furnsh(which('lsk/naif0011.tls'));
cspice_furnsh(which('pck/gm_de431.tpc'));
gmSun = cspice_bodvrd( 'Sun', 'GM', 1 ); % [km^3/s^2]
gmEarth = cspice_bodvrd( 'Earth', 'GM', 1 ); % [km^3/s^2]
gmMoon = cspice_bodvrd( 'Moon', 'GM', 1 ); % [km^3/s^2]

% define LU, TU, VU (mass unit is useless here), and mass parameter (mu)
%%% Earth-Moon system.
lengthUnit = 0.3844e9; % [m] % https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
timeUnit = sqrt( (lengthUnit/1e3)^3 / (gmEarth+gmMoon) ); % [s]
velocityUnit = lengthUnit / timeUnit; % [m/s]
mu = gmMoon / ( gmEarth + gmMoon );
disp('# In Earth-Moon system...');
%
%%% Sun-Earth/Moon system.
% lengthUnit = 149.60e9; % [m] % https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
% timeUnit = sqrt( (lengthUnit/1e3)^3 / (gmSun+gmEarth+gmMoon) ); % [s]
% velocityUnit = lengthUnit / timeUnit; % [m/s]
% mu = ( gmEarth + gmMoon ) / ( gmSun + gmEarth + gmMoon );
% disp('# In Sun-Earth/Moon system...');


%% generate Halo 3rd order approximation

location = 'L1';
phi = 0;

Az = 40000e3 / lengthUnit; % [m] --> [LU] % could only reach 40000e3 km for Earth-Moon now.

% segment information
segmentNumber = 52;

initialEpoches = linspace( 0, 20*pi, segmentNumber+1 );

% third order guess in CRTBP (segment number is determined here)
initialStatesThirdOrder = zeros(length(initialEpoches),6);
for ii = 1:length(initialEpoches)
    [initialStatesThirdOrder(ii,:), period] = VerticalLyapunovThirdOrder(location, phi, Az, mu, initialEpoches(ii));
end
% show 3rd order orbit
figure(93); 
% clf;
tt = linspace(initialEpoches(1),initialEpoches(end),1001);
X = zeros(length(tt),6);
for ii = 1:length(tt)
    X(ii,:) = VerticalLyapunovThirdOrder(location, phi, Az, mu, tt(ii));
end
plot3(X(:,1),X(:,2),X(:,3),'-'); hold on;
plot3(X(1,1),X(1,2),X(1,3),'b*'); hold on;
plot3(X(end,1),X(end,2),X(end,3),'ro'); hold on;
xlabel('x'); ylabel('y'); zlabel('z'); 
axis equal;


%% shooting in CRTBP
[initialStatesCRTBP,period] = VerticalLyapunovShooting(mu, initialStatesThirdOrder(1,:), period, 3, 1e-7);
odeOptions = odeset('AbsTol',1e-9,'RelTol',1e-9); 
[tempT,tempX] = ode113(@(t,X)DynamicRTBP(t,X,mu,0),...
    linspace( initialEpoches(1), initialEpoches(1)+period, 1001 ),...
    initialStatesCRTBP,...
    odeOptions); 
plot3(tempX(:,1),tempX(:,2),tempX(:,3));
% axis equal;
% interpolate resulted Halo in CRTBP
initialStates = interp1( tempT, tempX, mod(initialEpoches,period) );
%
% [another choice is to directly refine in the ephemeris model, but might need more segments]
% initialStates = initialStatesThirdOrder;

% show initial guess in crtbp
figure(93);
PlotInitialState(@(t,X)DynamicRTBP(t,X,mu,0), initialEpoches, initialStates);

%% multiple-shooting in moon-centered inertial ephemeris model
cspice_kclear;
spiceKernelList = {'lsk/naif0011.tls',  'spk/planets/de432s.bsp',  'pck/gm_de431.tpc', ...
    'EM_MoonCenteredRotation.fk',  'MoonCenteredInertial.fk'};
cspice_kclear;
for ii = 1:length(spiceKernelList)
    cspice_furnsh(which(spiceKernelList{ii}));
end

% nondimensionanl EMBR --> dimensional MCR
% after conversion, assumed in EMBR mice frame
%%% epoches
initialEphmerisEpoch = cspice_str2et('2017/11/13 00:00:00'); % [s] % the ephemeris epoch corresponding to initialEpoches(1)
initialEphmerisEpoches = initialEphmerisEpoch + (initialEpoches-initialEpoches(1)) * timeUnit; % [s] % the time interval matters
%%% states
initialEphemerisStatesECR = zeros(size(initialStates));
initialEphemerisStatesECR(:,1) = (initialStates(:,1)-(1-mu)) * lengthUnit / 1e3; % [km]
initialEphemerisStatesECR(:,2:3) = initialStates(:,2:3) * lengthUnit / 1e3; % [km]
initialEphemerisStatesECR(:,4:6) = initialStates(:,4:6) * lengthUnit / 1e3 / timeUnit; % [km/s]

% MCR --> MCI
TransMatrix = cspice_sxform('MCR','MCI',initialEphmerisEpoches);
initialEphemerisStatesECI = zeros(size(initialEphemerisStatesECR));
for ii = 1:length(initialEphmerisEpoches)
    initialEphemerisStatesECI(ii,:) = ( squeeze(TransMatrix(:,:,ii)) * initialEphemerisStatesECR(ii,:).' ).';
end

% dynamic equation in ephemeris model
fcn = @(t,X)DynamicEphemerisInertial(t,X,'Moon',{'Sun','Venus','Mercury','Earth','Moon','4','5','6','7','8'}, spiceKernelList);


% start shooting
positionTolerance = 0.01; % [km]
velocityTolerance = 0.0001; % [km/s]
odeOptions = odeset('AbsTol',1e-9,'RelTol',1e-9);
tic;
[correctedInitialEpoches, correctedInitialStates, exitflag] = MultipleShooting(fcn, initialEphmerisEpoches, initialEphemerisStatesECI, positionTolerance, velocityTolerance, odeOptions);
toc;

%% draw in dimensional results
figure(92); 
clf;
tic;
for ii = 1:length(correctedInitialEpoches)-1
    
    [plotT,plotX] = ode113(fcn,[correctedInitialEpoches(ii),correctedInitialEpoches(ii+1)],correctedInitialStates(ii,:),odeOptions);
    
    TransMatrix = cspice_sxform('MCI','MCR',plotT.');
    for jj = 1:length(plotT)
        plotX(jj,:) = ( squeeze(TransMatrix(:,:,jj)) * plotX(jj,:).' ).';
    end
    
    plot3(plotX(:,1),plotX(:,2),plotX(:,3),'-'); hold on;
    if ii == 1; plot3(plotX(1,1),plotX(1,2),plotX(1,3),'b*'); hold on; end
    if ii == length(correctedInitialEpoches)-1; plot3(plotX(end,1),plotX(end,2),plotX(end,3),'ro'); hold on; end
%     axis equal;
    xlabel('x [km]');
    ylabel('y [km]');
    zlabel('z [km]');
end
grid on;
toc;


%% generate orbit data

% in Moon Centered Inertial frame
[nominalTMCI, nominalXMCI] = PlotInitialState(fcn, correctedInitialEpoches, correctedInitialStates, odeOptions);
plot3(nominalXMCI(:,1),nominalXMCI(:,2),nominalXMCI(:,3),'b');

% in Moon centered Rotating frame
nominalXMCR = zeros(size(nominalXMCI));
TransMatrix = cspice_sxform('MCI','MCR',nominalTMCI.');
for jj = 1:length(nominalTMCI)
    nominalXMCR(jj,:) = ( squeeze(TransMatrix(:,:,jj)) * nominalXMCI(jj,:).' ).';
end
plot3(nominalXMCR(:,1),nominalXMCR(:,2),nominalXMCR(:,3),'r');

