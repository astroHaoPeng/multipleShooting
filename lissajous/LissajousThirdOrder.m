function X = LissajousThirdOrder(Position, phi, psi, Ay, Az, mu, t)
% Generate third order approximation of Lissajous orbit around L1 or L2 point. 
% (validated using multiple shooting method)
%
% Plan:
%   1. what about L3?
% 
% References:
%   Richardson, David L. 1980. “Analytic Construction of Periodic Orbits about the Collinear Points.” Celestial Mechanics 22 (3):241–253. https://doi.org/10.1007/BF01229511.
%   (need to change z-frequency and psi to use input values, other parameters are the same to halo's)
%
% see also: MultipleShooting, HaloThirdOrder
%
% created by PH at 2017-11-13:1113 based on HaloThirdOrder.m



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
c3 = c(3);
c4 = c(4);
lambda = roots([1 0 c2-2 0 -(c2-1)*(1+2*c2)]);
lambda = lambda( imag(lambda)==0 & lambda>0 );
k = 2*lambda/(lambda^2+1-c2);
nu = sqrt(c2);


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

Ax = Ay/k;

tau1 = lambda*t + phi;
tau2 = nu*t + psi;
switch Position % This is default to northern groups.
    case {'L1'}
        sigma = 1; % or -1, switch function
    case {'L2'}
        sigma = -1;
end

x = a21*Ax^2 + a22*Az^2 -Ax*cos(tau1) + (a23*Ax^2-a24*Az^2)*cos(2*tau1)...
    + (a31*Ax^3-a32*Ax*Az^2)*cos(3*tau1);
y = k*Ax*sin(tau1) + (b21*Ax^2-b22*Az^2)*sin(2*tau1)...
    + (b31*Ax^3-b32*Ax*Az^2)*sin(3*tau1);
z = sigma*(Az*cos(tau2)+d21*Ax*Az*(cos(2*tau2)-3)+(d32*Ax^2*Az-d31*Az^3)*cos(3*tau2));
vx = (Ax*sin(tau1) - 2*(a23*Ax^2-a24*Az^2)*sin(2*tau1) - 3*(a31*Ax^3-a32*Ax*Az^2)*sin(3*tau1))*lambda;
vy = (k*Ax*cos(tau1) + 2*(b21*Ax^2-b22*Az^2)*cos(2*tau1) + 3*(b31*Ax^3-b32*Ax*Az^2)*cos(3*tau1))*lambda;
vz = sigma*(-Az*sin(tau2)-2*d21*Ax*Az*sin(2*tau2)-3*(d32*Ax^2*Az-d31*Az^3)*sin(3*tau2))*nu;


X = [x,y,z,vx,vy,vz];
X = gamma.*X;
switch Position
    case {'L1'}
        X(1) = X(1) + 1 - mu - gamma;
    case {'L2'}
        X(1) = X(1) + 1 - mu + gamma;
end

end
