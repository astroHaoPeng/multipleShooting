addpath('halo');
addpath('rtbp');

%% define constants

%%% define LU, MU, TU
lengthUnit = 149598261e3; % [m] 
massUnit = 1.98855e30 + 5.97219e24 + 7.34767309e22; % [kg]
% timeUnit = 

%%% mass parameter
% Earth+Moon / Sun+Earth+Moon
mu = ( 5.97219e24 + 7.34767309e22 ) / ( 1.98855e30 + 5.97219e24 + 7.34767309e22 );


%% generate Halo in CRTBP

%
xL1 = LibrationPoint(mu, 'L1'); location = 'L1';
% xL2 = LibrationPoint(mu, 'L2'); location = 'L2';

% Halo orbit size Az
amplitudeZ = 120000e3 / lengthUnit; % [LU]
% 3rd approximation
initialPhase = pi;
X3 = HaloThirdOrder(amplitudeZ,'Az',location,mu,initialPhase);
% shooting method
tic;
[X0,~,~,exitflag] = HaloShooting( mu, X3, 5, 1e-9  ); % do not change [1,3,5] for current
toc;
disp(exitflag);

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