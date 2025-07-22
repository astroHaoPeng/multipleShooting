%Demonstration of generating Halo orbit in the Sun-Earth/Moon CRTBP
%
% `Demo_01_SunEarthMoon_Halo_CRTBP.m`: Generate Halo orbit in Sun-Earth CRTBP model using third-order approximation and simple shooting method.
%
% 1. 3rd order approximation of Halo
% 2. simple shooting, vary two of three nonzero components in 
%    initial contion[x0, 0, z0, 0, vy0, 0], but default we keep z0 fixed.
% 3. draw and show.
%
%Limits:
% see HaloShooting.m


%% add necessary paths
tmpExcludingExtension = {'.git', '.svn'}; % to exclude paths containing these extensions
tmpAllPaths = strsplit(genpath('../'), pathsep); % all paths
tmpAllPathsFiltered = strjoin(tmpAllPaths(~contains(tmpAllPaths, tmpExcludingExtension)), pathsep); % filtering the path
addpath(tmpAllPathsFiltered);

%% define constants

%%% define LU, MU, TU
gmSun = 1.3271e+11; % [km^3/s^2] % extracted from 'gm_de431.tpc', approximately
gmEarth = 3.9860e+05; % [km^3/s^2]
gmMoon = 4.9028e+03; % [km^3/s^2]
%%% Sun-Earth/Moon system.
lengthUnit = 149.60e9; % [m] % https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
timeUnit = sqrt( (lengthUnit/1e3)^3 / (gmSun+gmEarth+gmMoon) ); % [s]
velocityUnit = lengthUnit / timeUnit; % [m/s]
mu = ( gmEarth + gmMoon ) / ( gmSun + gmEarth + gmMoon );
mu = 0.00001
disp('# In Sun-Earth/Moon system...');


%% generate Halo in CRTBP

%
xL1 = LibrationPoint(mu, 'L1'); location = 'L1';
% xL2 = LibrationPoint(mu, 'L2'); location = 'L2';

% Halo orbit size Az
amplitudeZ = 220000e3 / lengthUnit; % [LU]
% 3rd approximation
initialPhase = pi;
[X3,initialPeriod] = HaloThirdOrder(amplitudeZ,'Az',location,mu,initialPhase);
% shooting method
tic;
[X0,~,~,exitflag] = HaloShooting( mu, X3, initialPeriod, 3, 1e-9  ); % do not change [1,3,5] for current
toc;
% disp(exitflag);

%% plot to test

figure(91);
clf

% 3rd order approximation
phi = initialPhase + linspace(0,0.98*2*pi,101); % [rad]
PV = zeros(length(phi),6); % [LU, VU]
for ii = 1:length(phi)
    PV(ii,:) = HaloThirdOrder(amplitudeZ,'Az',location,mu, phi(ii) );
end
p3rd = plot3(PV(:,1),PV(:,2),PV(:,3),'.:b'); hold on;

% shooting results
[~,PVHalo] = ode113(...
    @(f,X)DynamicRTBP(f,X,mu,0),...
    linspace(0,0.98*HaloPeriod(X0,mu),101),...
    X0,...
    odeset('AbsTol',1e-9,'RelTol',1e-9) );
pC = plot3(PVHalo(:,1),PVHalo(:,2),PVHalo(:,3),'.-r'); hold on;
pCS = plot3(PVHalo(1,1),PVHalo(1,2),PVHalo(1,3),'*b', 'markersize',15); hold on;
pCE = plot3(PVHalo(end,1),PVHalo(end,2),PVHalo(end,3),'or', 'markersize',15); hold on;

% Libration point
if strcmpi(location,'L1')
    pL = plot3(xL1, 0, 0, 'kx', 'markersize',15); hold on;
elseif strcmpi(location,'L2')
    pL = plot3(xL2, 0, 0, 'kx', 'markersize',15); hold on;
end
axis equal;

% labels
if strcmpi(location,'L1')
    legend([pL, p3rd, pC, pCS, pCE], 'L_1', '3rd-order approx.', 'shooting in CRTBP', 'start', 'end', 'Location','northeast');
elseif strcmpi(location,'L2')
    legend([pL, p3rd, pC, pCS, pCE], 'L_2', '3rd-order approx.', 'shooting in CRTBP', 'start', 'end', 'Location','northeast');
end

set(gca,'FontSize',12);
grid on;
%%