function Xdot = DynamicRTBP(f,X,mu,e)
% Restricted Three-Body Problem 动力学方程
%
%   以真近点角为独立变量；CRTBP 时即为时间，ERTBP 时为角度
%
%   调用格式：
%       Xdot = RTBPOde(f,X,mu,e,number_vectors,delta_f)
%
%   Input:
%       f,X,mu,e:   
%       number_vectors: 向量化计算时 X 中包含的点数
%       delta_f:        向量化的各段与输入的时刻 f 之间的时间差：向量化要求各段的时间相等，f+delta_f 得到各段的真实时间
%
%   Output:
%       Xdot:       与输入的 X 相同维数的向量
%
%
%   注意：勿轻易修改，小的改动会导致小的数值误差，先前计算的轨道库的精度就会受到影响
%
% 	created by PH at 2014-04-21:1842: 从 HaloOde 修改而来，以后都调用这个
% 	last modified by PH at 2014-08-12:1247 完善注释；增加对向量输入的支持（未验证）
%   PH at 2014-11-06:1436 需要对多输入进行验证
%
%   created by PH at 2017-11-12:1920 
%       recreated based on RTBPOde.m for testing multiple shooting project
%       remove support for vectorized input

%% equation of dynamic motion
% 由 symbolic.m 推导来的ODE
x = X(1);
y = X(2);
z = X(3);
dx = X(4);
dy = X(5);
dz = X(6);

r1cubic = ((mu + x).^2 + y.^2 + z.^2).^(3/2);
r2cubic = ((mu + x - 1).^2 + y.^2 + z.^2).^(3/2);
if e == 0
    kappa = 1;
else
    kappa = (e.*cos(f) + 1);
end

Xdot = [
    dx;
    dy;
    dz;
     2*dy + (x + ((2*mu + 2*x).*(mu - 1))./(2*r1cubic) - (mu.*(2*mu + 2*x - 2))./(2*r2cubic))./kappa;
    -2*dx + (y - (mu*y)./r2cubic + (y.*(mu - 1))./r1cubic)./kappa;
    -z + (z - (mu*z)./r2cubic + (z.*(mu - 1))./r1cubic)./kappa];

%% first order vairation equation
if length(X) > 6
    r1power5 = ((mu + x)^2 + y^2 + z^2)^(5/2);
    r2power5 = ((mu + x - 1)^2 + y^2 + z^2)^(5/2);
    two_x = 2 * x;
    two_mu = 2 * mu;
    three_mu = 3 * mu;

    for ii = 1:6
        phix = X(ii*6+1);
        phiy = X(ii*6+2);
        phiz = X(ii*6+3);
        phivx = X(ii*6+4);
        phivy = X(ii*6+5);
        phivz = X(ii*6+6);
    
        Xdot = [...
            Xdot;
            phivx;
            phivy;
            phivz;
            2*phivy + (phix*((mu-1)/r1cubic - mu/r2cubic + (three_mu*(two_mu + two_x - 2)^2)/(4*r2power5) - (3*(two_mu + two_x)^2*(mu-1))/(4*r1power5) + 1))/kappa + (phiy*((three_mu*y*(two_mu + two_x - 2))/(2*r2power5) - (3*y*(two_mu + two_x)*(mu-1))/(2*r1power5)))/kappa + (phiz*((three_mu*z*(two_mu + two_x - 2))/(2*r2power5) - (3*z*(two_mu + two_x)*(mu-1))/(2*r1power5)))/kappa;
            (phiy*((mu-1)/r1cubic - mu/r2cubic - (3*y^2*(mu-1))/r1power5 + (three_mu*y^2)/r2power5 + 1))/kappa - 2*phivx - (phiz*((3*y*z*(mu-1))/r1power5 - (three_mu*y*z)/r2power5))/kappa + (phix*((three_mu*y*(two_mu + two_x - 2))/(2*r2power5) - (3*y*(two_mu + two_x)*(mu-1))/(2*r1power5)))/kappa;
            phiz*(((mu-1)/r1cubic - mu/r2cubic - (3*z^2*(mu-1))/r1power5 + (three_mu*z^2)/r2power5 + 1)/kappa - 1) - (phiy*((3*y*z*(mu-1))/r1power5 - (three_mu*y*z)/r2power5))/kappa + (phix*((three_mu*z*(two_mu + two_x - 2))/(2*r2power5) - (3*z*(two_mu + two_x)*(mu-1))/(2*r1power5)))/kappa ];
    end
end

end