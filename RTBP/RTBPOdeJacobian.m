% 给出 RTBP 的动力学方程的 Jacobian 矩阵
%   在计算状态转移矩阵，或者当前点的线性稳定性时使用
%   与当前点的线性化等价
%
%   created by PH at 2014-06-19:1258 由
%       E:\Codes\MATLAB\[____已经完成的程序____]\20130930 重新编写Halo轨道程序 的
%       symbolic2_modified 推导而来

function A = RTBPOdeJacobian(f,X,mu,e)

x = X(1);
y = X(2);
z = X(3);
dx = X(4);
dy = X(5);
dz = X(6);

r2cubic = ((mu + x - 1)^2 + y^2 + z^2)^(3/2);
r1cubic = ((mu + x)^2 + y^2 + z^2)^(3/2);
kappa = (e*cos(f) + 1);
r2power5 = ((mu + x - 1)^2 + y^2 + z^2)^(5/2);;
r1power5 = ((mu + x)^2 + y^2 + z^2)^(5/2);

A = [
    [0,0,0,  1, 0, 0]
    [0,0,0,  0, 1, 0]
    [0,0,0,  0, 0, 1]
    [ ((mu - 1)/r1cubic - mu/r2cubic + (3*mu*(2*mu + 2*x - 2)^2)/(4*r2power5) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*r1power5) + 1)/kappa,      ((3*mu*y*(2*mu + 2*x - 2))/(2*r2power5) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*r1power5))/kappa,          ((3*mu*z*(2*mu + 2*x - 2))/(2*r2power5) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*r1power5))/kappa,     0, 2, 0]
    [                                     ((3*mu*y*(2*mu + 2*x - 2))/(2*r2power5) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*r1power5))/kappa,    ((mu - 1)/r1cubic - mu/r2cubic - (3*y^2*(mu - 1))/r1power5 + (3*mu*y^2)/r2power5 + 1)/kappa,                                           -((3*y*z*(mu - 1))/r1power5 - (3*mu*y*z)/r2power5)/kappa,    -2, 0, 0]
    [                                     ((3*mu*z*(2*mu + 2*x - 2))/(2*r2power5) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*r1power5))/kappa,                                       -((3*y*z*(mu - 1))/r1power5 - (3*mu*y*z)/r2power5)/kappa,    ((mu - 1)/r1cubic - mu/r2cubic - (3*z^2*(mu - 1))/r1power5 + (3*mu*z^2)/r2power5 + 1)/kappa - 1,     0, 0, 0]
    ];

end
