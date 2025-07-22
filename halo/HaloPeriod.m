function period = HaloPeriod(X0, mu, odeOptions)
%Calculate period of Halo orbit defined by X0.
% X0 should be a point on Halo. But it can be any point on the Halo. 
% So `HaloShooting` must be first used be generate a perfect Halo orbit.
%
%   1. Find the previous cross with x-z plane;
%   2. Propagate until the next cross.
%   3. Double the propagated time, which is the period.
%
%   period = HaloPeriod(X0, mu)
%   period = HaloPeriod(X0, mu, odeOptions)
%
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
%private function to stop propagation when corss with x-z plane, where y == 0.
value = X(2);
isterminal = 1;
direction = 0;
end