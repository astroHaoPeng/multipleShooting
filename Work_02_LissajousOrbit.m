%% define constants

%%% define LU, MU, TU
lengthUnit = 149598261e3; % [m] 
massUnit = 1.98855e30 + 5.97219e24 + 7.34767309e22; % [kg]
timeUnit = 365.25*86400/2/pi; % [s]
velocityUnit = lengthUnit / timeUnit; % [m/s]


%%% mass parameter
% Earth+Moon / Sun+Earth+Moon
mu = ( 5.97219e24 + 7.34767309e22 ) / ( 1.98855e30 + 5.97219e24 + 7.34767309e22 );

%% generate first order Lissajous
Position = 'L1';
phi = pi;
psi = pi/2;
Ay = 300000e3 / lengthUnit; % [km]
Az = 150000e3 / lengthUnit; % [km]

tt = linspace(0,2.5*2*pi,101);
for ii = 1:length(tt)
    X(ii,:) = LissajousFirstOrder(Position, phi, psi, Ay, Az, mu, tt(ii));
end

figure(93);
clf
plot3(X(:,1),X(:,2),X(:,3),'.-'); hold on;
xlabel('x');
ylabel('y');
zlabel('z');
% axis equal;

% %% multiple shooting in CRTBP
initialEpoches = linspace( 0, 2.5*2*pi, 21 );
initialStates = zeros(length(initialEpoches),6);
for ii = 1:length(initialEpoches)
    initialStates(ii,:) = LissajousFirstOrder(Position, phi, psi, Ay, Az, mu, initialEpoches(ii));
end
% PlotInitialState(@(t,X)DynamicRTBP(t,X,mu,0), initialEpoches, initialStates, mu);
tic;
[correctedInitialEpoches, correctedInitialStates, exitflag] = MultipleShooting(@(t,X)DynamicRTBP(t,X,mu,0), initialEpoches, initialStates);
toc;
% PlotInitialState(correctedInitialEpoches, correctedInitialStates, mu);




%%
function X = LissajousFirstOrder(Position, phi, psi, Ay, Az, mu, t)
% 计算L1位置 L2位置
switch Position
    case {'L1'}
        gamma = roots([1 -(3-mu) 3-2*mu -mu  2*mu -mu]);
        gamma = gamma(imag(gamma)==0);
    case {'L2'}
        gamma = roots([1  (3-mu) 3-2*mu -mu -2*mu -mu]);
        gamma = gamma(imag(gamma)==0);
end

% 重新规一划
Ay = Ay/gamma;
Az = Az/gamma;

% 一阶近似
switch Position
    case {'L1'}
        c = @(n)(        mu + (-1)^n*(1-mu)*gamma^(n+1)/(1-gamma)^(n+1) ) / gamma^3;
    case {'L2'}
        c = @(n)( (-1)^n*mu + (-1)^n*(1-mu)*gamma^(n+1)/(1+gamma)^(n+1) ) / gamma^3;
end
c2 = c(2);
lambda = roots([1 0 c2-2 0 -(c2-1)*(1+2*c2)]);
lambda = lambda( imag(lambda)==0 & lambda>0 );
k = 2*lambda/(lambda^2+1-c2);
nu = sqrt(c2);


x = -Ay*cos(lambda*t+phi);
y = k*Ay*sin(lambda*t+phi);
z = Az*sin(nu*t+psi);
vx = Ay*sin(lambda*t+phi)*lambda;
vy = k*Ay*cos(lambda*t+phi)*lambda;
vz = Az*cos(nu*t+psi)*nu;

X = [x,y,z,vx,vy,vz];
X = gamma.*X;
switch Position
    case {'L1'}
        X(1) = X(1) + 1 - mu - gamma;
    case {'L2'}
        X(1) = X(1) + 1 - mu + gamma;
end

end
