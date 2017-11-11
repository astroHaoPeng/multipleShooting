function Phi = HaloPhi(tspan, X0, PhiColumn, mu, e, Tol)
%% 计算对Halo三阶近似进行微分校正需要的Phi矩阵
% Phi = HaloPhi(tspan, X0, PhiColumn, mu, e, Tol)
% last modified by PH at 2013-10-12:1046 提高积分精度，之前积分精度小于微分修正精度，不合理
% last modified by PH at 2013-10-21:2156 引入Tol来控制精度，之前的调用都需要修改
% last modified by PH at 2013-10-21:2156 修改了对PhiColumn，使之可以用parfor并行计算


%% 输入检测
%%% 默认偏心率e=0
if nargin == 4
    e = 0;
end
if nargin == 5
    Tol.RelTol = 1e-12;
    Tol.AbsTol = 1e-12;
    Tol.MaxStep = 1e-3;
end
if nargin == 6
    if ~isfield(Tol,'RelTol'); Tol.RelTol = 1e-12; end;
    if ~isfield(Tol,'AbsTol'); Tol.AbsTol = 1e-12; end;
    if ~isfield(Tol,'MaxStep'); Tol.MaxStep = 1e-3; end;
end
%%% 如果输入tspan为一个值，则区间为[0,tspan]，兼容之前的程序
if length(tspan) == 1
    tspan = [0, tspan];
end

%% 计算需要的状态转移矩阵
Phi = zeros(6);
finiteDifference = 1e-6;
options = odeset('RelTol',Tol.RelTol,'AbsTol',Tol.AbsTol);
options = odeset(options, 'MaxStep', Tol.MaxStep);
for ii = 1:6
    if any( ii == PhiColumn )
        Y0 = [zeros(6,1); reshape(X0,6,1)];
        Y0(ii) = finiteDifference;
        [~, Y] = ode113(@(t,Y)trans(t,Y,mu,e), tspan, Y0, options);
        Phi(:,ii) = reshape(Y(end,1:6), 6, 1);
    end
end

%% 输出Phi
Phi = Phi / finiteDifference;
end

function dY = trans(f,Y,mu,e)
xx = Y(1);
yy = Y(2);
zz = Y(3);
dxx = Y(4);
dyy = Y(5);
dzz = Y(6);
x = Y(7);
y = Y(8);
z = Y(9);
% dx = Y(10);
% dy = Y(11);
% dz = Y(12);

% dY = [
%     dxx
%     dyy
%     dzz
%     2*dyy + (xx*((mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1))/(e*cos(f) + 1) + (yy*((3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))))/(e*cos(f) + 1) + (zz*((3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))))/(e*cos(f) + 1)
%     (yy*((mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1))/(e*cos(f) + 1) - 2*dxx - (zz*((3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) - (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2)))/(e*cos(f) + 1) + (xx*((3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))))/(e*cos(f) + 1)
%     (zz*((mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1))/(e*cos(f) + 1) - (yy*((3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) - (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2)))/(e*cos(f) + 1) + (xx*((3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))))/(e*cos(f) + 1)
%     reshape(HaloOde(f,Y(7:12),mu,e), 6, 1)];

dY = [dxx
    dyy
    dzz
    2*dyy + (xx*((mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2 + z^2)^(5/2)) + 1))/(e*cos(f) + 1) + (yy*((3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))))/(e*cos(f) + 1) + (zz*((3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))))/(e*cos(f) + 1)
    (yy*((mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1))/(e*cos(f) + 1) - 2*dxx - (zz*((3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) - (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2)))/(e*cos(f) + 1) + (xx*((3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))))/(e*cos(f) + 1)
    zz*(((mu - 1)/((mu + x)^2 + y^2 + z^2)^(3/2) - mu/((mu + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) + (3*mu*z^2)/((mu + x - 1)^2 + y^2 + z^2)^(5/2) + 1)/(e*cos(f) + 1) - 1) - (yy*((3*y*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(5/2) - (3*mu*y*z)/((mu + x - 1)^2 + y^2 + z^2)^(5/2)))/(e*cos(f) + 1) + (xx*((3*mu*z*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2 + z^2)^(5/2))))/(e*cos(f) + 1)
    reshape(DynamicRTBP(f,Y(7:12),mu,e), 6, 1)];

end
