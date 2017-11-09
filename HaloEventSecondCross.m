function [value,isterminal,direction] = HaloEventSecondCross(f,X,X0)
% [value,isterminal,direction] = HaloEventSecondCross(f,X,Position)
% 第二次与x-z平面相交
%
% last modified by PH at 2013-10-12:1046
% last modified by PH at 2013-10-29:1033 加入了对L2点的支持，Position
% last modified by PH at 2013-11-27:1958
%   修改了对 direction 的赋值逻辑，与Position无关，第二次相次，应该与初始值的方向相同即可。
%   要求输入初始值 X0
% last modified by PH at 2013-11-28:2043 更新了根据 X0 对位置的判断，可以兼容 Halo 的左右两侧
% last modified by PH at 2013-11-28:2043 更新为基于 HaloEventFirsrCross.m，方便修改


% 输入检测
if ischar(X0)
    error('We have changed: use X0 to replace Position value!');
end

%
[value,isterminal,direction] = HaloEventFirstCross(f,X,X0);
direction = -direction;

% % 以 y 为事件值
% value = X(2);
% 
% % 终止
% isterminal = 1;
% 
% % 根据 X0 判断对应的 Halo 轨道的位置
% % 以及是位于 Halo 轨道与 x-z 平面的交点是在 Left 还是 Right
% if X0(1)<1.05 && X0(3)>0 && X0(5)>0
%     direction = 1; % 第二次相交
% elseif X0(1)<1.05 && X0(3)<0 && X0(5)<0
%     direction = -1;
% elseif X0(1)>0.95 && X0(3)<0 && X0(5)>0
%     direction = 1;
% elseif X0(1)>0.95 && X0(3)>0 && X0(5)<0
%     direction = -1;
% end

end