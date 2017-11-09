function Period = HaloPeriod(t0,X0,mu,Tol)
% 计算给定初始条件 t0, X0 的 Halo 轨道的周期
%   要求输入的 X0 必须在 x-z 平面上，即 X0([2,4,6]) == 0
%
%   last modified by PH at 2014-06-26:1437
%   last modified by PH at 2014-07-08:1058 暂时归入函数库，还需要改进，判断最终停止积分时的精度

% 输入检测
% warning('In library, but needs more modification...')

% if any( X0([2,4,6]) ~= 0 )
%     error('X0(2), X0(4), X0(6) must be zeros');
% end
if ~exist('Tol','var') || isempty(Tol)
    Tol.AbsTol = 1e-7;
    Tol.RelTol = 1e-7;
    warning('Default Tolerance used...');
end

% 找到周期 Period
OdeOptions = odeset('AbsTol',Tol.AbsTol, 'RelTol',Tol.RelTol);
OdeOptions = odeset(OdeOptions, 'Event', @(t,X)HaloEventSecondCross(t,X,X0));
[tempT, ~] = ode45(@(t,X)DynamicRTBP(t,X,mu,0), [t0, 2*pi], X0, OdeOptions);
Period = tempT(end) - t0;

end