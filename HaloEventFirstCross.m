function [value,isterminal,direction] = HaloEventFirstCross(f,X,X0)
% [value,isterminal,direction] = HaloEventFirstCross(f,X,X0)
% 第一次与x-z平面相交
% 
% last modified by PH at 2013-10-12:1046
% last modified by PH at 2013-11-27:1711 修改了注释部分，加入调用方法
% last modified by PH at 2013-11-27:1958
% 修改了对 direction 的赋值逻辑，与Position无关，第二次相次，应该与初始值的方向相同即可。
% 要求输入初始值 X0
% last modified by PH at 2013-11-28:2033 更新了根据 X0 对位置的判断，可以兼容 Halo 的左右两侧


% 输入检测
if ischar(X0)
    error('We have changed: use X0 to replace Position value!');
end

% 以 y 为事件值
value = X(2);

% 终止
isterminal = 1;

% 根据 X0 判断对应的 Halo 轨道的位置
% 以及是位于 Halo 轨道与 x-z 平面的交点是在 Left 还是 Right
if X0(1)<1.05 && X0(3)>0 && X0(5)>0
    direction = -1; % 第一次相交，则 dy/dt 方向应该和初始值的 dy/dt 相反
elseif X0(1)<1.05 && X0(3)<0 && X0(5)<0
    direction = 1;
elseif X0(1)>0.95 && X0(3)<0 && X0(5)>0
    direction = -1;
elseif X0(1)>0.95 && X0(3)>0 && X0(5)<0
    direction = 1;
end

end