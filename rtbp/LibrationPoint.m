function [x, y] = LibrationPoint(mu, Position)
%% [x, y] = LibrationPoint(mu, Position)
% 根据输入mu计算L1或者L2点的坐标，规一旋转坐标系下表示
%
% last modified by PH at 2013-10-12:1046
% last modified by PH at 2013-10-12:1046
% 加入了L3,L4,L5的支持，并通过max和min添加了对mu=0和mu=1的支持，输出更改为两位坐标

if mu < 0 || mu > 1
    error('Mass ratio mu must be 0<=mu<=1');
end

% Parker & Anderson, Low-energy Lunar Trajectory Design, 2013
switch Position
    case {'L1','l1'}
        gamma = roots([1 -(3-mu) 3-2*mu -mu  2*mu -mu]);
        gamma = gamma(imag(gamma)==0);
        x = max( 1 - mu - gamma );
        y = 0;
    case {'L2','l2'}
        gamma = roots([1  (3-mu) 3-2*mu -mu -2*mu -mu]);
        gamma = gamma(imag(gamma)==0);
        x = max( 1 - mu + gamma );
        y = 0;
    case {'L3','l3'}
        gamma = roots([1 (2+mu) (1+2*mu) -(1-mu) -2*(1-mu) -(1-mu)]);
        gamma = gamma(imag(gamma)==0);
        x = min( - mu - gamma );
        y = 0;
    case {'L4','l4'}
        x = 0.5 - mu;
        y = sqrt(3)/2;
    case {'L5','l5'}
        x = 0.5 - mu;
        y = - sqrt(3)/2;
    otherwise
        error('Please provide correct libration point Position.');
end
end