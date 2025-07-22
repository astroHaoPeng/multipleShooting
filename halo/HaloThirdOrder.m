function [X0, period] = HaloThirdOrder(A, Direction, Position, mu, phi, t)
%function X0 = HaloThirdOrder(A, Direction, Position, mu, phi)
% Generates a third-order analytical approximate solution for the Halo orbit
%   A           Amplitude
%   Dirction    Dirction amplitude direction, either 'Ax' or 'Az'
%   Position    L1 or L2
%   mu
%   phi         Initial phase on the orbit, 0 for the left, pi for the right (where see from above)

% Modified by PH
% Modified Date: 20130930
% last modified by PH at 2014-01-06:1414 加入了参数 phi 的支持
%   当 phi=0 时给出左边起点的初始值
%   当 phi=pi 时给出右边起点的初始值，分别对应于 Left 和 Right Groups 的需求

if nargin == 0
    mu = 3.040357143e-6; % earth-moon and sun system
    L_e2s = 1.495978714e8; % distance from the earth to the sun
    A = 145000/L_e2s;
    Direction = 'Az';
    Position = 'L1';
elseif nargin == 3 % When phi is not given, the default is the starting point on the left, which is consistent with the previous call
    phi = 0; 
end

%% Calculate L1 position L2 position
switch Position
    case {'L1'}
        gamma = roots([1 -(3-mu) 3-2*mu -mu  2*mu -mu]);
        gamma = gamma(imag(gamma)==0);
    case {'L2'}
        gamma = roots([1  (3-mu) 3-2*mu -mu -2*mu -mu]);
        gamma = gamma(imag(gamma)==0);
end

%% Re-regularization
A = A/gamma;

%% Calculate the coefficients of the third-order approximate solution, according to the reference:
% David L. Richardson, "Analytic Construction of Periodic Orbits about the Collinear Points”, Celestial Mechanics, vol. 22, Oct. 1980, pp. 241–253.
switch Position
    case {'L1'}
        c = @(n)(        mu + (-1)^n*(1-mu)*gamma^(n+1)/(1-gamma)^(n+1) ) / gamma^3;
    case {'L2'}
        c = @(n)( (-1)^n*mu + (-1)^n*(1-mu)*gamma^(n+1)/(1+gamma)^(n+1) ) / gamma^3;
end
c2 = c(2);
c3 = c(3);
c4 = c(4);
lambda = roots([1 0 c2-2 0 -(c2-1)*(1+2*c2)]);
lambda = lambda(logical((imag(lambda)==0).*(lambda>0)));
k = 2*lambda/(lambda^2+1-c2);
period = 2*pi/lambda; % <--------- This is the 3rd-order approximated period.

d1 = 3*lambda^2/k*(k*(6*lambda^2-1)-2*lambda);
d2 = 8*lambda^2/k*(k*(11*lambda^2-1)-2*lambda);

a21 = 3*c3*(k^2-2)/4/(1+2*c2);
a22 = 3*c3/4/(1+2*c2);
a23 = -3*c3*lambda/4/k/d1*(3*k^3*lambda-6*k*(k-lambda)+4);
a24 = -3*c3*lambda/4/k/d1*(2+3*k*lambda);
b21 = -3*c3*lambda/2/d1*(3*k*lambda-4);
b22 = 3*c3*lambda/d1;
d21 = -c3/2/lambda^2;

a31 = -9*lambda/4/d2*(4*c3*(k*a23-b21)+k*c4*(4+k^2))...
    + (9*lambda^2+1-c2)/2/d2*(3*c3*(2*a23-k*b21)+c4*(2+3*k^2));
a32 = -1/d2*(9*lambda/4*(4*c3*(k*a24-b22)+k*c4)...
    +3/2*(9*lambda^2+1-c2)*(c3*(k*b22+d21-2*a24)-c4));
b31 = 3/8/d2*(8*lambda*(3*c3*(k*b21-2*a23)-c4*(2+3*k^2))...
    +(9*lambda^2+1+2*c2)*(4*c3*(k*a23-b21)+k*c4*(4+k^2)));
b32 = 1/d2*(9*lambda*(c3*(k*b22+d21-2*a24)-c4)...
    +3/8*(9*lambda^2+1+2*c2)*(4*c3*(k*a24-b22)+k*c4));
d31 = 3/64/lambda^2*(4*c3*a24+c4);
d32 = 3/64/lambda^2*(4*c3*(a23-d21)+c4*(4+k^2));

a1 = -3/2*c3*(2*a21+a23+5*d21)-3/8*c4*(12-k^2);
a2 = 3/2*c3*(a24-2*a22)+9/8*c4;
s1 = 1/(2*lambda*(lambda*(1+k^2)-2*k))*(3/2*c3*(2*a21*(k^2-2)-a23*(k^2+2)-2*k*b21)...
    -3/8*c4*(3*k^4-8*k^2+8));
s2 = 1/(2*lambda*(lambda*(1+k^2)-2*k))*(3/2*c3*(2*a22*(k^2-2)+a24*(k^2+2)+2*k*b22+5*d21)...
    +3/8*c4*(12-k^2));
l1 = a1 + 2*lambda^2*s1;
l2 = a2 + 2*lambda^2*s2;
Delta = lambda^2 - c2;

%% The third-order approximate solution for generating Halo orbits depending on whether `Ax` or `Az` is given
switch Direction
    case {'Xmax'}
        Ax = A-1;
        Az = sqrt((-Delta-l1*Ax^2)/l2);
    case {'Ax'}
        Ax = A;
        Az = sqrt((-Delta-l1*Ax^2)/l2);
    case {'Az'}
        Az = A;
        Ax = sqrt((-Delta-l2*Az^2)/l1);
end
if ~exist('t','var')
    t = 0;
end
tau1 = lambda*t + phi;
switch Position % The default is the north halo oribt, whose part above the x-z plane is longer.
    case {'L1'}
        sigma = 1; % or -1, switch function
    case {'L2'}
        sigma = -1;
end

x = a21*Ax^2 + a22*Az^2 -Ax*cos(tau1) + (a23*Ax^2-a24*Az^2)*cos(2*tau1)...
    + (a31*Ax^3-a32*Ax*Az^2)*cos(3*tau1);
y = k*Ax*sin(tau1) + (b21*Ax^2-b22*Az^2)*sin(2*tau1)...
    + (b31*Ax^3-b32*Ax*Az^2)*sin(3*tau1);
z = sigma*(Az*cos(tau1)+d21*Ax*Az*(cos(2*tau1)-3)+(d32*Ax^2*Az-d31*Az^3)*cos(3*tau1));
Vx = (Ax*sin(tau1) - 2*(a23*Ax^2-a24*Az^2)*sin(2*tau1) - 3*(a31*Ax^3-a32*Ax*Az^2)*sin(3*tau1))*lambda;
Vy = (k*Ax*cos(tau1) + 2*(b21*Ax^2-b22*Az^2)*cos(2*tau1) + 3*(b31*Ax^3-b32*Ax*Az^2)*cos(3*tau1))*lambda;
Vz = sigma*(-Az*sin(tau1)-2*d21*Ax*Az*sin(2*tau1)-3*(d32*Ax^2*Az-d31*Az^3)*sin(3*tau1))*lambda;

%% Output the initial value of the third-order approximate Halo orbit and restore it to the CR3BP coordinate system
X0 = [x; y; z; Vx; Vy; Vz];
X0 = gamma.*X0;
switch Position
    case {'L1'}
        X0(1) = X0(1) + 1 - mu - gamma;
    case {'L2'}
        X0(1) = X0(1) + 1 - mu + gamma;
end

end