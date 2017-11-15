addpath('ephemeris');
addpath('rtbp');

clear;

parpool;

%% define constants

%%% define LU, MU, TU
lengthUnit = 149598261e3; % [m] 
massUnit = 1.98855e30 + 5.97219e24 + 7.34767309e22; % [kg]
timeUnit = 365.25*86400/2/pi; % [s]
velocityUnit = lengthUnit / timeUnit; % [m/s]


%%% mass parameter
% Earth+Moon / Sun+Earth+Moon
mu = ( 5.97219e24 + 7.34767309e22 ) / ( 1.98855e30 + 5.97219e24 + 7.34767309e22 );

%% generate Lissajous initial guess for shooting
Position = 'L1';
phi = pi;
psi = pi/2;
Ay = 300000e3 / lengthUnit; % [km]
Az = 150000e3 / lengthUnit; % [km]
% Ay = 700000e3 / lengthUnit; % [km]
% Az = 900000e3 / lengthUnit; % [km]


% third order guess in CRTBP (segment number is determined here)
segmentNumber = 20;
initialEpoches = linspace( 0, 2.5*2*pi, segmentNumber+1 );
initialStates = zeros(length(initialEpoches),6);
for ii = 1:length(initialEpoches)
    initialStates(ii,:) = LissajousThirdOrder(Position, phi, psi, Ay, Az, mu, initialEpoches(ii));
end

% % shooting in CRTBP
% tic;
% [initialEpoches, initialStates, exitflag] = MultipleShooting(@(t,X)DynamicRTBP(t,X,mu,0), initialEpoches0, initialStates0);
% toc;

% show 3rd order orbit
figure(93);
clf
tt = linspace(initialEpoches(1),initialEpoches(end),1001);
X = zeros(length(tt),6);
for ii = 1:length(tt)
    X(ii,:) = LissajousThirdOrder(Position, phi, psi, Ay, Az, mu, tt(ii));
end
plot3(X(:,1),X(:,2),X(:,3),'.-'); hold on;
plot3(X(1,1),X(1,2),X(1,3),'b*'); hold on;
plot3(X(end,1),X(end,2),X(end,3),'ro'); hold on;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

% show initial guess in crtbp
% figure(93);
% PlotInitialState(@(t,X)DynamicRTBP(t,X,mu,0), initialEpoches, initialStates);


%% convert to earth

cspice_kclear;
miceRootFolder = '/Users/GroupMacBai/Codes/Tools/SPICE/generic_kernels/';
cspice_furnsh([miceRootFolder 'lsk/naif0011.tls']);
cspice_furnsh([miceRootFolder 'spk/planets/de432s.bsp']);
cspice_furnsh([miceRootFolder 'pck/gm_de431.tpc']);
cspice_furnsh('EarthCenteredRotation.fk');
cspice_furnsh('EarthCenteredInertial.fk');

% nondimensionanl --> dimensional
% after conversion, assumed in EMBR mice frame
%%% epoches
initialEphmerisEpoch = cspice_str2et('2017/11/13 00:00:00'); % [s] % the ephemeris epoch corresponding to initialEpoches(1)
initialEphmerisEpoches = initialEphmerisEpoch + (initialEpoches-initialEpoches(1)) * timeUnit; % [s] % the time interval matters
%%% states
initialEphemerisStatesECR = zeros(size(initialStates));
initialEphemerisStatesECR(:,1) = (initialStates(:,1)-(1-mu)) * lengthUnit / 1e3; % [km]
initialEphemerisStatesECR(:,2:3) = initialStates(:,2:3) * lengthUnit / 1e3; % [km]
initialEphemerisStatesECR(:,4:6) = initialStates(:,4:6) * lengthUnit / 1e3 / timeUnit; % [km/s]

% EMBR --> EMBI
TransMatrix = cspice_sxform('ECR','ECI',initialEphmerisEpoches);
initialEphemerisStatesECI = zeros(size(initialEphemerisStatesECR));
for ii = 1:length(initialEphmerisEpoches)
    initialEphemerisStatesECI(ii,:) = ( squeeze(TransMatrix(:,:,ii)) * initialEphemerisStatesECR(ii,:).' ).';
end

% draw orbit in ECI frame
% fcn = @(t,X)DynamicEphemerisInertial(t,X,'Earth',{'Sun','Earth','Moon'});
fcn = @(t,X)DynamicEphemerisInertial(t,X,'Earth',{'Sun','Venus','Mercury','Earth','Moon','4','5','6','7','8'});
PlotInitialState(fcn, initialEphmerisEpoches, initialEphemerisStatesECI, odeset('AbsTol',1e-9,'RelTol',1e-9) );

%% shooting
positionTolerance = 0.01; % [km]
velocityTolerance = 0.0001; % [km/s]
odeOptions = odeset('AbsTol',1e-6,'RelTol',1e-6);
tic;
[correctedInitialEpoches, correctedInitialStates, exitflag] = MultipleShooting(fcn, initialEphmerisEpoches, initialEphemerisStatesECI, positionTolerance, velocityTolerance, odeOptions);
toc;

%%
figure(92); 
clf;
tic;
for ii = 1:length(correctedInitialEpoches)-1
    
    [plotT,plotX] = ode113(fcn,[correctedInitialEpoches(ii),correctedInitialEpoches(ii+1)],correctedInitialStates(ii,:));
    
    TransMatrix = cspice_sxform('ECI','ECR',plotT.');
    for jj = 1:length(plotT)
        plotX(jj,:) = ( squeeze(TransMatrix(:,:,jj)) * plotX(jj,:).' ).';
    end
    
%     subplot(121)
    plot3(plotX(:,1),plotX(:,2),plotX(:,3),'-'); hold on;
    if ii == 1; plot3(plotX(1,1),plotX(1,2),plotX(1,3),'b*'); hold on; end
    if ii == length(correctedInitialEpoches)-1; plot3(plotX(end,1),plotX(end,2),plotX(end,3),'ro'); hold on; end
    axis equal;
    xlabel('x');
    ylabel('y');
    
%     subplot(122)
%     plot3(plotX(:,3),plotX(:,4),plotX(:,5),'-'); hold on;
%     plot3(plotX(1,3),plotX(1,4),plotX(1,5),'b*'); hold on;
%     plot3(plotX(end,3),plotX(end,4),plotX(end,5),'ro'); hold on;
%     xlabel('vx');
%     ylabel('vy');
end
toc;


%% show results
% figure(93);
% PlotInitialState(@(t,X)DynamicRTBP(t,X,mu,0), correctedInitialEpoches, correctedInitialStates);





