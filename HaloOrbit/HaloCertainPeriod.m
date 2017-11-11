function [tempHalo,ii,direction] = HaloCertainPeriod(T, OrbitDatabase, Direction, Position, mu, e, Tol)
%% 采用二分法，在给定的轨道库OrbitDatabase中，找到周期为T的Halo轨道
% tempHalo = HaloEllipticInitial(T, OrbitDatabase, Direction, Position, mu, e)
%
% Last modified by PH at 2013-10-15:1836 改变了名字，原来是HaloEllipticInitial
% Last modified by PH at 2013-10-28:2118 修改了判断和二分，使可以兼容单调增的情况
%                                        ！！还不支持二次型最小值的情况，需要注意
% Last modified by PH at 2013-10-28:2118 加入了精度结构Tol控制周期精度
% Last modified by PH at 2014-07-08:1114 修正了输出，直接以 XLeft 为输出，与计算 err 时相符合


% 输入检测
if nargin == 6
    Tol.PeriodTol = 1e-5;
    Tol.RelTol = 1e-13;
    Tol.AbsTol = 1e-13;
    Tol.CorrectionTol = 1e-11;
end  
if ~isfield(Tol,'PeriodTol') % 周期确定精度
    Tol.PeriodTol = 1e-5;
end

%%% 找到在轨道库中的位置
N = length(OrbitDatabase(:,1));
for ii = 1:N
    if T < OrbitDatabase(ii, 7) && T > OrbitDatabase(ii+1, 7)
        direction = -1;
        break; % 递减时退出
    end
    if T > OrbitDatabase(ii, 7) && T < OrbitDatabase(ii+1, 7)
        direction = 1;
        break; % 递增时退出
    end
end

% 用二分法来找足够精度的周期为T的轨道
XLeft = OrbitDatabase(ii, 1:10);
XRight = OrbitDatabase(ii+1, 1:10);
iter = 1; % 二分计数器
err = abs(XLeft(7)-T);
while err > Tol.PeriodTol % 以左边的周期为目标周期的近似
    timer1 = tic;
    
    % 直接二分
% 	XMiddle = HaloContinuation((XLeft(10)+XRight(10))/2, Direction, Position,...
%         XLeft(1:6), abs(-XLeft(10)+XRight(10))/2, mu, e, Tol); 
		
	% 按比例二分
    AzNew = XLeft - (XLeft(7)-T)/(XLeft(7)-XRight(7)) * (XLeft-XRight);
    AzNew = AzNew(10);
    AzStep = abs( AzNew - XLeft(10) );
    XMiddle = HaloContinuation(AzNew, Direction, Position, XLeft(1:6), AzStep, mu, e, Tol);

    switch direction
        case 1 % 递增
            if T < XMiddle(7)
                XRight = XMiddle;
            elseif T > XMiddle(7)
                XLeft = XMiddle;
            else
                break;
            end
        case -1 % 递减
            if T < XMiddle(7)
                XLeft = XMiddle;
            elseif T > XMiddle(7)
                XRight = XMiddle;
            else
                break;
            end
    end
    %%% 输出二分次数
    err = abs(XLeft(7)-T);
    fprintf('二分法第%3d 次，用时%6.3f，方向是%2d，\n  TLeft=%10.9f, TRight=%10.9f, err=%10.9f\n',...
        iter, toc(timer1), direction, XLeft(7), XRight(7), err);
    iter = iter + 1;
    
end

%% 输出周期为T的Halo轨道的初始值
% tempHalo = [X0, T, Ax, Ay, Az]
% % 选择距离T较近的一端为输出值
% if (XLeft(7)-T) > (T-XRight(7))
%     tempHalo = XRight;
% else
%     tempHalo = XLeft;
% end
% % 选择线性加权值为输出值
% tempHalo = XLeft - (XLeft(7)-T)/(XLeft(7)-XRight(7)) * (XLeft-XRight);
% 直接以 XLeft 为输出，因为在计算 err 时是以 XLeft(7) 为准
tempHalo = XLeft;

end
