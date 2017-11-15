function period = HaloPeriod(X0, mu, odeOptions)
%Calculate period of Halo orbit defined by X0.
%period = HaloPeriod(X0, mu)
%period = HaloPeriod(X0, mu, odeOptions)
%
% X0 should be a point on Halo. But it can be any point on the Halo. 
% So HaloShooting should be first used be generate a perfect Halo orbit. 
%
%
% 计算给定初始条件 t0, X0 的 Halo 轨道的周期
%   X0 可以是 Halo 轨道上的任意点。
%
% see also: HaloShooting

if nargin < 4
    odeOptions = odeset('RelTol',1e-9, 'AbsTol',1e-9);
end

odeOptions = odeset(odeOptions, 'Events',@(t,X)PrivateFirstCross(t,X,X0));

%% backwardly find the first cross with x-y plane
if abs(X0(2)) > 1e-9
    % find start point
    [tempT, tempX] = ode113(@(t,X)DynamicRTBP(t,X,mu,0), [0, -2*pi], X0, odeOptions);
    t0 = tempT(end);
else
    t0 = 0;
end

%% fowardly find the first cross
% find end point
[tempT, tempX] = ode113(@(t,X)DynamicRTBP(t,X,mu,0), [0, 2*pi], X0, odeOptions);
t1 = tempT(end);

%% Then period is twice of the duration
period = 2 * (t1 - t0);

end



function [value,isterminal,direction] = PrivateFirstCross(f,X,X0)
%private function to stop propagation when corss with x-z plan, where y == 0.
value = X(2);
isterminal = 1;
direction = 0;
end