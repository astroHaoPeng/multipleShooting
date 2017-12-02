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
    jacobian = zeros(36,1);
    for ii = 1:6
        jacobian((ii-1)*6+1:(ii-1)*6+3) = X(ii*6+4:ii*6+6);
    end
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
%     P = cspice_spkezr( targ, et, 'J2000', 'None', inertialFrameCenter ); % [km, km/s] % state of gravity body
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
    % first order vairation equation
    if length(X) > 6
        px = P(1);
        py = P(2);
        pz = P(3);
        Rpower3 = ((px - x)^2 + (py - y)^2 + (pz - z)^2)^(3/2);
        Rpower5 = ((px - x)^2 + (py - y)^2 + (pz - z)^2)^(5/2);
        threeGM = 3*GM(ii);
        
        for jj = 1:6
            phix = X(jj*6+1);
            phiy = X(jj*6+2);
            phiz = X(jj*6+3);
            %phivx = X(jj*6+4); % not used
            %phivy = X(jj*6+5);
            %phivz = X(jj*6+6);
            
            jacobian((jj-1)*6+4:(jj-1)*6+6) = jacobian((jj-1)*6+4:(jj-1)*6+6)...
                + [...
                    (threeGM*phiy*(py - y)*(px - x))/(Rpower5) - phix*(GM(ii)/Rpower3 - (threeGM*(px - x)*(px - x))/(Rpower5)) + (threeGM*phiz*(pz - z)*(px - x))/(Rpower5)
                    (threeGM*phix*(px - x)*(py - y))/(Rpower5) - phiy*(GM(ii)/Rpower3 - (threeGM*(py - y)*(py - y))/(Rpower5)) + (threeGM*phiz*(pz - z)*(py - y))/(Rpower5)
                    (threeGM*phix*(px - x)*(pz - z))/(Rpower5) - phiz*(GM(ii)/Rpower3 - (threeGM*(pz - z)*(pz - z))/(Rpower5)) + (threeGM*phiy*(py - y)*(pz - z))/(Rpower5)
                    ];
        end
    end
end



%% output
Xdot = [X(4:6); a];
if length(X) == 42
    Xdot = [Xdot; jacobian];
end

end