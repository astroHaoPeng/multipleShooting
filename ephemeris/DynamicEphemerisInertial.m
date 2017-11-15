function Xdot = DynamicEphemerisInertial(et,X,inertialFrameCenter,gravityBodies)
%% 地月星历动力学模型，地月质心，惯性系 J2000 下进行积分
%
%   地月质心 J2000 惯性坐标系
%   
%   Input:
%       t:          [TU] 归一化时间 
%       X:          [LengthUnit, LengthUnit/TimeUnit] 地月质心 J2000 惯性坐标系（非脉动）
%       LengthUnit: [km] 输入 X 的归一化距离单位
%       TimeUnit:   [s] 输入 t 的归一化时间单位
%       varargin:   需要考虑的 JPL 星历天体
%   
%   Ouput:
%       dX:         输入 X 的单位相同
%
%   created by PH at 2014-07-21:1135
%   last modified by PH at 2014-08-12:1338 完善注释
%   last modified by Ph at 2014-08-12:1404 复制到函数库（未测试）
% warning('untested!');

if cspice_ktotal('all') ~= 5
    cspice_kclear;
    miceRootFolder = '/Users/GroupMacBai/Codes/Tools/SPICE/generic_kernels/';
    cspice_furnsh([miceRootFolder 'lsk/naif0011.tls']);
    cspice_furnsh([miceRootFolder 'spk/planets/de432s.bsp']);
    cspice_furnsh([miceRootFolder 'pck/gm_de431.tpc']);
    cspice_furnsh('EarthCenteredRotation.fk');
    cspice_furnsh('EarthCenteredInertial.fk');
end


if size(X,2)>1
    error('X should be a row vector.');
end

if length(X)>6
    jacobian = zeros(36,1);
    for ii = 1:6
        jacobian((ii-1)*6+1:(ii-1)*6+3) = X(ii*6+4:ii*6+6);
    end
end

x = X(1);
y = X(2);
z = X(3);

a = 0; % [km/s^2]
for ii = 1:length(gravityBodies)
    
    % body parameters
    targ = gravityBodies{ii};
    GM = cspice_bodvrd( targ, 'GM', 1 ); % [km^3/s^2]
    P = cspice_spkezr( targ, et, 'J2000', 'None', inertialFrameCenter ); % [km, km/s] % state of gravity body 
    
    % gravity acceloration
    if strcmpi( targ, inertialFrameCenter ) 
        % center body
        a = a - GM / norm( X(1:3) - P(1:3) )^3 * ( X(1:3) - P(1:3) ); % [km/s] 
    else
        % perturbation body
        a = a - GM / norm( X(1:3) - P(1:3) )^3 * ( X(1:3) - P(1:3) )...
            - GM / norm( P(1:3) )^3 * P(1:3); % [km/s^2]
    end
    
    % first order vairation equation
    if length(X) > 6
        px = P(1);
        py = P(2);
        pz = P(3);
        Rpower3 = ((px - x)^2 + (py - y)^2 + (pz - z)^2)^(3/2);
        Rpower5 = ((px - x)^2 + (py - y)^2 + (pz - z)^2)^(5/2);
        threeGM = 3*GM;
        
        for ii = 1:6
            phix = X(ii*6+1);
            phiy = X(ii*6+2);
            phiz = X(ii*6+3);
            %phivx = X(ii*6+4);
            %phivy = X(ii*6+5);
            %phivz = X(ii*6+6);
            
            jacobian((ii-1)*6+4:(ii-1)*6+6) = jacobian((ii-1)*6+4:(ii-1)*6+6)...
                + [...
                    (threeGM*phiy*(py - y)*(px - x))/(Rpower5) - phix*(GM/Rpower3 - (threeGM*(px - x)*(px - x))/(Rpower5)) + (threeGM*phiz*(pz - z)*(px - x))/(Rpower5)
                    (threeGM*phix*(px - x)*(py - y))/(Rpower5) - phiy*(GM/Rpower3 - (threeGM*(py - y)*(py - y))/(Rpower5)) + (threeGM*phiz*(pz - z)*(py - y))/(Rpower5)
                    (threeGM*phix*(px - x)*(pz - z))/(Rpower5) - phiz*(GM/Rpower3 - (threeGM*(pz - z)*(pz - z))/(Rpower5)) + (threeGM*phiy*(py - y)*(pz - z))/(Rpower5)
                    ];
        end
    end
end



%% output
Xdot = [X(4:6); a];
if length(X)>6
    Xdot = [Xdot; jacobian];
end

end