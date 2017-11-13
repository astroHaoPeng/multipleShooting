function dX = RTBP_EarthMoon_Ephemeris(t,X,LengthUnit,TimeUnit,varargin)

%% 星历动力学模型
%
%   dX = RTBP_EarthMoon_Ephemeris(t,X,LengthUnit,TimeUnit,varargin)
%
%   地月 RTBP 模型
%
%   Input:
%       t:          归一化时间 [TU]
%       X:          state [LengthUnit, LengthUnit/TimeUnit]
%       LengthUnit: [km] 输入 X 的归一化距离单位
%       TimeUnit:   [s] 输入 t 的归一化时间单位
%       varargin:   需要考虑的 JPL 星历天体
%   
%   Ouput:
%       dX:         输入 X 的单位相同
%
%   created by PH at 2014/07/01 17:02
%   last modified by PH at 2014-07-20:1430 完全修改了程序，按照归一化时间来计算
%   last modified by PH at 2014-08-12:1338 完善注释
%

%%
error('This model is wrong!!! Use EphemerisEarthMoon.m as instead.')

et = t * TimeUnit; % [s] 星历时间

% 地月旋转系下的摄动加速度
a = 0; % [km/s^2]
ref = 'EMBR';
obs = 'EARTH MOON BARYCENTER';
for ii = 1:length(varargin)
    targ = varargin{ii};
    % 摄动行星的位置
    state_P = cspice_spkezr( targ, et, ref, 'None', obs ); % [km, km/s] % Planet
    GM_P = cspice_bodvrd( targ, 'GM', 1 ); % [km^3/s^2]
    a = a - GM_P / norm( X(1:3)*LengthUnit - state_P(1:3) )^3 * ( X(1:3)*LengthUnit - state_P(1:3) )...
        + GM_P / norm( state_P(1:3) )^3 * state_P(1:3); % [km/s^2]
end

% 地月产生的加速度
ref = 'EMBR';
obs = 'EARTH MOON BARYCENTER';
state_E = cspice_spkezr( 'Earth', et, ref, 'None', obs ); % [km, km/s] % Earth
state_M = cspice_spkezr( 'Moon',  et, ref, 'None', obs ); % [km, km/s] % Moon
GM_E = cspice_bodvrd( 'Earth', 'GM', 1 ); % [km^3/s^2]
GM_M = cspice_bodvrd( 'Moon',  'GM', 1 ); % [km^3/s^2]
r_E_SC = X(1:3)*LengthUnit - state_E(1:3);
r_M_SC = X(1:3)*LengthUnit - state_M(1:3);
a = a - GM_E / norm( r_E_SC )^3 * r_E_SC...
    - GM_M / norm( r_M_SC )^3 * r_M_SC;

% 地月系瞬时角速度，归一化
ref = 'EMBI';
obs = 'EARTH MOON BARYCENTER';
state_E_inertial = cspice_spkezr( 'Earth', et, ref, 'None', obs ); % [km, km/s] % Earth
omega = cross( state_E_inertial(1:3), state_E_inertial(4:6) ) ...
    / norm( state_E_inertial(1:3) )^2 ...
    * TimeUnit; % [rad, s^-1]

% 归一化加速度
a = a / LengthUnit * TimeUnit^2;

% 科氏力和旋转坐标系
a = a + -2 * cross(omega,X(4:6)) - cross( omega, cross(omega,X(1:3)) );

dX = [X(4:6); a];

% % disp(et);
% disp(norm(omega));

end