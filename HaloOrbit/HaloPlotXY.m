function [T,X] = HaloPlotXY(X0, mu, e, tspan, RELTOL, ABSTOL, MaxStep)
%% 绘制Halo轨道XY平面投影图
% HaloPlotXY(X0, mu, e, tspan, RELTOL, ABSTOL)
%
% last modified by PH at 2013-10-17:2106
% last modified by PH at 2014-04-21:1838 功能已经被 HaloPlot 取代，此处不再更新


%% 输入检测和默认参数
%%% 绘图精度
if nargin <= 4
    RELTOL = 1e-11;
    ABSTOL = 1e-11;
end
if nargin <= 6
    MaxStep = [];
end
odeOptions = odeset('RelTol',RELTOL,'AbsTol',ABSTOL);
odeOptions = odeset(odeOptions, 'MaxStep',MaxStep);
%%% 输入检测
if nargin <= 3 % 此时不输入tspan，可能人为确定
    odeOptions = odeset('RelTol',RELTOL,'AbsTol',ABSTOL, 'Events',@HaloEventSecondCross);
    tspan = [0,10];
end
if nargin == 2
    e = 0;
end

%% 绘图
[T X] = ode45(@(t,X)HaloOde(t,X,mu,e), tspan, X0, odeOptions);

plot(X(:,1),X(:,2));
grid on; hold on; axis equal;
plot(X(1,1),X(1,2),'b*');
plot(X(end,1),X(end,2),'ro');
title('x-y');

end