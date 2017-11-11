function X1 = HaloPoint(t0,X0,t1,Period,mu,e,Tol)
% 根据 t0 时的位置 X0 计算 t1 时的位置 X1
%
%   X1 = HaloPoint(t0,X0,t1,Period,mu,e,Tol)
%
% created by PH at 2014-01-09:1535
% last modified by PH at 2014-06-26:1424 添加了默认精度，和 warning
% last modified by PH at 2014-07-08:1054 归入函数库

% 输入检测
if ~exist('Tol','var') || isempty(Tol)
    Tol.AbsTol = 1e-7;
    Tol.RelTol = 1e-7;
    warning('Default Tolerance used...');
end

% 被 Period 约化后的轨道上的点
% 这样可以保证计算精度
t1hat = mod(t1,Period);

if t1hat == t0
    X1 = reshape(X0,1,6);
else
    OdeOptions = odeset('AbsTol',Tol.AbsTol, 'RelTol',Tol.RelTol);
    [~,tempX] = ode45(@(t,X)HaloOde(t,X,mu,e), [t0,t1hat], X0, OdeOptions);
    X1 = tempX(end,:);
end

end
