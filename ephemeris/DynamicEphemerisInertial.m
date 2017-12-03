function Xdot = DynamicEphemerisInertial(et,X,inertialFrameCenter,gravityBodies,spiceKernels)
%% ephemeris model in inertial frame parallel with J2000
%   centered at 'inertialFrameCenter' (string)
%   gravatiational attraction from all 'gravityBodies' (string cell) are considered
%   NAIF SPICE toolkit is used with kernels in 'spiceKernels' (string cell)
%   What not considered: solar radiation, ...
%
%   Input:
%       et: ephemeris time, usually in [s]
%       X:  ephemeris state, [x,y,z,vx,vy,vz], usually in [km] and [km/s], which is consistent with SPICE
%       inertialFrameCenter: center of the inertial frame
%       gravityBodies: as the variable name says, should be covered by 'spiceKernels'
%       spiceKernels: {'kernel_01',...} kernels used by SPICE
%   
%   Ouput:
%       Xdot: acceleration, usually in [km/s^2]
%
%   created by PH at 2017-12-01 17:44 UTC-4 @Rutgers:  
%   updated by PH at 2017-12-03 16:05 UTC-4 @Rutgers: update STM calculation, faster now

%% input check and initialization

% clear cache and load all kernels at the first call
persistent flagKernelLoaded;
if isempty(flagKernelLoaded)
    cspice_kclear;
    % load all kernels from input
    for ii = 1:length(spiceKernels)
        cspice_furnsh(which(spiceKernels{ii}));
    end
    flagKernelLoaded = 1;
end

% creat jacobian if input X is 42
if length(X)==42
    jacobianLowerLeft = zeros(3,3);
elseif length(X)~=6
    error('Input X should be 6x1');
end

%%
x = X(1);
y = X(2);
z = X(3);

a = 0; % [km/s^2]

persistent GM;

for ii = 1:length(gravityBodies)
    
    %--------------------------------------
    % extract body parameters: GM, position
    targ = gravityBodies{ii};
    if length(GM)<ii || isempty(GM(ii)) % using persistent variables, just extract all GM for only one, which can greatly speedup
        GM(ii) = cspice_bodvrd( targ, 'GM', 1 ); % [km^3/s^2]
    end
    P = cspice_spkpos( targ, et, 'J2000', 'None', inertialFrameCenter ); % [km, km/s] % state of gravity body
    % P has been tested, cannot be further speedup by using persistent variables

    %---------------------
    % gravity acceloration
    if strcmpi( targ, inertialFrameCenter ) 
        % center body
        a = a - GM(ii) / norm( X(1:3) - P(1:3) )^3 * ( X(1:3) - P(1:3) ); % [km/s^2] 
    else
        % perturbation body
        a = a - GM(ii) / norm( X(1:3) - P(1:3) )^3 * ( X(1:3) - P(1:3) )...
            - GM(ii) / norm( P(1:3) )^3 * P(1:3); % [km/s^2]
    end
    
    %-------------------------------
    % Jacobian matrix, or first order vairation equation, or Phi, or STM (state transition matrix)
    % The same for center body and perturbation body ! Because only position of satellite is varied.
    % It should be in the form:
    %   0   0   0   1   0   0
    %   0   0   0   0   1   0
    %   0   0   0   0   0   1
    %   Uxx Uxy Uxz 0   0   0
    %   Uyx Uyy Uyz 0   0   0
    %   Uzx Uzy Uzz 0   0   0
    % so, only the lower left corner is accumulated through different gravitational bodies
    if length(X) == 42
        px = P(1); % planet x
        py = P(2);
        pz = P(3);
        Rpower3 = ((px - x)^2 + (py - y)^2 + (pz - z)^2)^(3/2); % term to speedup
        Rpower5 = ((px - x)^2 + (py - y)^2 + (pz - z)^2)^(5/2); % term to speedup
        threeGM = 3*GM(ii); % term to speedup
        
        % Uxx--Uzz in lower left corner
        jacobianLowerLeft(1,1) = jacobianLowerLeft(1,1) + (threeGM*(px - x)*(px - x))/(Rpower5) - GM(ii)/Rpower3;
        jacobianLowerLeft(1,2) = jacobianLowerLeft(1,2) + (threeGM*(py - y)*(px - x))/(Rpower5);
        jacobianLowerLeft(1,3) = jacobianLowerLeft(1,3) + (threeGM*(pz - z)*(px - x))/(Rpower5);
        jacobianLowerLeft(2,1) = jacobianLowerLeft(1,2);
        jacobianLowerLeft(2,2) = jacobianLowerLeft(2,2) + (threeGM*(py - y)*(py - y))/(Rpower5) - GM(ii)/Rpower3;
        jacobianLowerLeft(2,3) = jacobianLowerLeft(2,3) + (threeGM*(pz - z)*(py - y))/(Rpower5);
        jacobianLowerLeft(3,1) = jacobianLowerLeft(1,3);
        jacobianLowerLeft(3,2) = jacobianLowerLeft(2,3);
        jacobianLowerLeft(3,3) = jacobianLowerLeft(3,3) + (threeGM*(pz - z)*(pz - z))/(Rpower5) - GM(ii)/Rpower3;
    end
    
end



%% output
if length(X) == 6
    % only acceloration
    Xdot = [X(4:6); a];
    
elseif length(X) == 42
    % acceloration and state transition matrix
    jacobian = [zeros(3),          eye(3);
                jacobianLowerLeft, zeros(3)];
    stateTransitionMatrix = jacobian * reshape(X(7:42),6,6);
    Xdot = [X(4:6); a; reshape(stateTransitionMatrix,36,1)];
    
end

end