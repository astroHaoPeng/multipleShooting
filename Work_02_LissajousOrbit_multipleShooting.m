addpath('lissajous');
addpath('rtbp');

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
% smaller
Ay = 300000e3 / lengthUnit; % [km]
Az = 150000e3 / lengthUnit; % [km]
% larger
Ay = 700000e3 / lengthUnit; % [km]
Az = 900000e3 / lengthUnit; % [km]

% generate 3rd order orbit
tt = linspace(0,2.5*2*pi,1001);
X = zeros(length(tt),6);
for ii = 1:length(tt)
    X(ii,:) = LissajousThirdOrder(Position, phi, psi, Ay, Az, mu, tt(ii));
end

% show 3rd order orbit
figure(93);
clf
plot3(X(:,1),X(:,2),X(:,3),'.-'); hold on;
plot3(X(1,1),X(1,2),X(1,3),'b*'); hold on;
plot3(X(end,1),X(end,2),X(end,3),'ro'); hold on;
xlabel('x');
ylabel('y');
zlabel('z');
% axis equal;

% initial guess in CRTBP (segment number is determined here)
segmentNumber = 15;
initialEpoches = linspace( 0, 2.5*2*pi, segmentNumber+1 );
initialStates = zeros(length(initialEpoches),6);
for ii = 1:length(initialEpoches)
    initialStates(ii,:) = LissajousThirdOrder(Position, phi, psi, Ay, Az, mu, initialEpoches(ii));
end

% show initial guess in crtbp
% figure(93);
% PlotInitialState(@(t,X)DynamicRTBP(t,X,mu,0), initialEpoches, initialStates);

%% shooting
tic;
[correctedInitialEpoches, correctedInitialStates, exitflag] = MultipleShooting(@(t,X)DynamicRTBP(t,X,mu,0), initialEpoches, initialStates);
toc;

%% show results
figure(93);
PlotInitialState(@(t,X)DynamicRTBP(t,X,mu,0), correctedInitialEpoches, correctedInitialStates);





