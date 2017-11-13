function [tempHalo,iter] = HaloContinuation(A,Direction,Position, X0,step, mu, e, Tol)
%% 把Halo轨道从初始条件X0延拓到振幅为Direction方向的A，步长为step
% tempHalo = HaloContinuation(A,Direction,Position, X0,step, mu, e, RELTOL, ABSTOL)
%
% last modified by PH at 2013-10-12:1046
% last modified by PH at 2013-10-23:1935 修正沿Ax延拓为沿Xmax延拓，即距m2天体的最远距离
% last modified by PH at 2013-10-28:2208 加入了精度控制Tol
% last modified by PH at 2013-10-30:0925 输出迭代次数 iter


%% 测试条件
if nargin == 0
    A = 8.3e-4;
    Direction = 'Az';
    Position = 'L1';
    A0 = 7.4e-4;
    step = 2e-6;
    mu = 3.040357143e-6; % earth-moon and sun system
    e = 0;
    X0 = HaloThirdOrder(A0, Direction, Position, mu);
    Tol.RelTol = 1e-11;
    Tol.AbsTol = 1e-11;
end

if nargin == 7
    Tol.RelTol = 1e-11;
    Tol.AbsTol = 1e-11;
    Tol.CorrectionTol = 1e-9;
end

%% 初始化输出
% X0 = [];
T = [];
Ax = [];
Ay = [];
Az = [];

%% 延拓
switch Direction
    case {'Az'}
        while abs( A - X0(3) ) > step
            if A > X0(3)
                X0(3) = X0(3) + step;
            else
                X0(3) = X0(3) - step;
            end
            [X0,iter] = HaloShooting(X0, [1,5], mu, e, Position, Tol); % 保证Az不变，即z不变
            if isempty(X0)
                tempHalo = [];
                return
            end
        end
        if A ~= X0(3)
            X0(3) = A;
            [X0,iter] = HaloShooting(X0, [1,5], mu, e, Position, Tol); % 保证Az不变，即z不变
            if isempty(X0)
                tempHalo = [];
                return
            end
        end
    case {'Xmax'} % 此时 A == Xmax
        tempA = 1 - mu - A;
        while abs( tempA - X0(1) ) > step
            if tempA > X0(1)
                X0(1) = X0(1) + step;
            else
                X0(1) = X0(1) - step;
            end
            [X0,iter] = HaloShooting(X0, [3,5], mu, e, Position, Tol); % 保证tempA，即Xmax
            if isempty(X0)
                tempHalo = [];
                return
            end
        end
        if tempA ~= X0(1)
            X0(1) = tempA;
            [X0,iter] = HaloShooting(X0, [3,5], mu, e, Position, Tol); % 保证tempA，即Xmax
            if isempty(X0)
                tempHalo = [];
                return
            end
        end
end



%% 生成输出变量
options = odeset('RelTol',Tol.RelTol, 'AbsTol',Tol.AbsTol, 'Events',@(f,X)HaloEventSecondCross(f,X,X0));
[T X] = ode45(@(t,X)HaloOde(t,X,mu,e), [0,3*pi], X0, options);
% Halo轨道周期
T = T(end);
% 计算各方向振幅，以 规一划单位 给出
switch Direction
    case {'Az'}
        Xmax = 1 - mu - X(1,1);
        Ay = max(X(:,2));
        Az = A;
    case {'Xmax'}
        Xmax = A;
        Ay = max(X(:,2));
        Az = X(1,3); % 初值的z振幅
end

%%% 输出为一个数组
tempHalo = [reshape(X0,1,6), T, Xmax, Ay, Az]; 

% HaloPlot(X0, mu, e); grid on; axis equal; hold on; % 画出轨道进行测试

end

