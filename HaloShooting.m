function [X0, iter] = HaloShooting(X0, PhiColumn, mu, e, Position, Tol)
% [X0, iter] = HaloShooting(X0, PhiColumn, mu, e, Position, Tol)
%
% inputs:
%   X0: initial guess. Usually from 3rd order approximation.
%   PhiColumn: usually [1,3,5]. This is index of columns that will be modifed during the shootting
%   mu:
%   e:
%   Position: 'L1' or 'L2'
%   Tol: tolerances used for shooting, must include all 4 fields
%       Tol.RelTol = 1e-13; (default)
%       Tol.AbsTol = 1e-13; (default)
%       Tol.CorrectionTol = 1e-9; (default)
%       Tol.MaxIteration = 30; (default)
% outputs:
%   X0: modifed state
%   iter: iteration number before exit. If iter == Tol.MaxIteration, it is very possible shooting has failed.
%
%% Halo轨道微分修正
%
% Last modified by PH at 2013-10-09:1547
% last modified by PH at 2013-10-12:1046 加入了最大迭代次数的限制，当超过最大迭代次数时，输出空数组
% last modified by PH at 2013-10-23:2225 加入一维线性搜索
% last modified by PH at 2013-10-28:2208 加入了精度控制Tol
% last modified by PH at 2013-10-29:1047 加入了Position控制，以适应L2的延拓
% last modified by PH at 2013-10-29:1500 增加了输出iter


%% 输入检测
%%% 默认精度和迭代次数
if nargin <= 5
    Tol.RelTol = 1e-13;
    Tol.AbsTol = 1e-13;
    Tol.CorrectionTol = 1e-9;
    Tol.MaxIteration = 30;
end
if ~isfield(Tol,'MaxIteration')
    Tol.MaxIteration = 10; % 默认迭代次数100次
end
if ~isfield(Tol,'CorrectionTol')
    Tol.CorrectionTol = 1e-9; % 默认修正精度1e-9
end

%%% 默认偏心率
if nargin <=3
    e = 0;
end

%% differential correction for the Halo orbit
odeOptions = odeset('RelTol',Tol.RelTol, 'AbsTol',Tol.AbsTol, 'Events',@(f,X)HaloEventFirstCross(f,X,X0));
[T X] = ode113(@(t,X)DynamicRTBP(t,X,mu,e), [0,pi], X0, odeOptions);
THalf = T(end);
X1 = X(end,:);
% plot(X(:,1), X(:,2),'r');grid on;axis equal;hold on;
% plot3(X(:,1), X(:,2), X(:,3), 'r');grid on;axis equal;hold on;

iter = 1;
while abs(X1(4))>Tol.CorrectionTol || abs(X1(6))>Tol.CorrectionTol
    Phi = HaloPhi(THalf, X0, PhiColumn, mu, e);
    temp = DynamicRTBP(THalf, X1, mu, e);
    L = [Phi(2,PhiColumn(1)), Phi(2,PhiColumn(2)), temp(2);
        Phi(4,PhiColumn(1)), Phi(4,PhiColumn(2)), temp(4);
        Phi(6,PhiColumn(1)), Phi(6,PhiColumn(2)), temp(6)];
    b = [X1(2);
        X1(4);
        X1(6)];
    delta = -L\b; % * (0.6+0.4*rand)
    s = fminsearch(@(s)tempfun(s,X0,delta,PhiColumn,mu,e,odeOptions), 1); % 一维线性搜索最优步长s，收敛性更好
%     s = 1; % 固定步长为 1，速度要快一点
    X0(PhiColumn(1)) = X0(PhiColumn(1)) + delta(1) * s;
    X0(PhiColumn(2)) = X0(PhiColumn(2)) + delta(2) * s;
    [T X] = ode113(@(t,X)DynamicRTBP(t,X,mu,e), [0,2*pi], X0, odeOptions);
    THalf = T(end);
    X1 = X(end,:);
%         plot(X(:,1), X(:,2));grid on;axis equal;hold on;pause(0.2);
    %     plot3(X(:,1), X(:,2), X(:,3));grid on;axis equal;hold on;
    if iter > Tol.MaxIteration
        X0 = [];
        return;
        warning('shooting exit with too many iterations');
    else
        iter = iter + 1;
    end
end


%  plot(X(:,1), X(:,2), 'g'); grid on;axis equal;hold on;

end


function Err = tempfun(s,X0,delta,PhiColumn,mu,e,odeOptions)
% 一维线性搜索的目标函数

X0(PhiColumn(1)) = X0(PhiColumn(1)) + delta(1) * s;
X0(PhiColumn(2)) = X0(PhiColumn(2)) + delta(2) * s;

[~, X] = ode113(@(t,X)DynamicRTBP(t,X,mu,e), [0,pi], X0, odeOptions);
Err = norm(X(end,2:2:6)); % 最终误差
end