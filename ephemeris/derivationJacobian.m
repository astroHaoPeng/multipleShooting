%% Derive the Jacobian matrix in a Barycenter Inertial Frame, meaning that there 
% is no massive body at the berycenter. Otherwise, the `aCenterPlanet` will be 
% derived differently.

%% Initialize symbols
syms GM
syms x y z vx vy vz
syms px py pz
syms phix phiy phiz phivx phivy phivz

% spacecraft
X = [x,y,z].'; % S/C position
v = [vx, vy, vz].'; % S/C velocity

% gravitational planet
P = [px,py,pz].'; % planet position


%% Acceleration due to a third planet

% (a of the S/C)  minus  (a of barycenter)
aOtherPlanet = - GM / sum( (X-P).^2 )^(3/2) * (X-P) - GM / sum( P.^2 )^(3/2) * P

% find jacobian of dX over all states
dX = [v; aOtherPlanet];
jacobianOtherPlanet = jacobian(dX, [x,y,z,vx,vy,vz]);

% get an expression that can be copy-and-pasted to write the ODE code.
phiX = [phix,phiy,phiz,phivx,phivy,phivz].';
simplify( jacobian(dX, [x,y,z,vx,vy,vz]) * phiX );


%% Similarly, the acceleration due to the centeral bodies

% (a of the S/C)  minus  (a of barycenter)
aCenterPlanet = - GM / sum( (X-P).^2 )^(3/2) * (X-P) - GM / sum( P.^2 )^(3/2) * P

% find jacobian of dX over all states
dX = [v; aCenterPlanet];
jacobianCenterPlanet = jacobian(dX, [x,y,z,vx,vy,vz]);

% get an expression that can be copy-and-pasted to write the ODE code.
phiX = [phix,phiy,phiz,phivx,phivy,phivz].';
simplify( jacobian(dX, [x,y,z,vx,vy,vz]) * phiX );


%% Verify that the ephemeris dynamics are the same for the third planets and centeral planets,
% as they are essentially the same situation.
% if there is a body in the center of the frame, the derivation will be different.

simplify( jacobianOtherPlanet - jacobianCenterPlanet )
