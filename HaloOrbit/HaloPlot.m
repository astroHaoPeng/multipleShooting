function [TAll,XAll, TPlot,XPlot] = HaloPlot(X0, mu, e, tspan, PlotType, Tol, MaxPlotPoint)
%[TAll,XAll, TPlot,XPlot] = HaloPlot(X0, mu, e, tspan, PlotType, Tol)
% 绘制 Ellipitc Halo 轨道，在脉动旋转系下
%
% 输入：
%  PlotType: 'x-y' 'y-z' 'x-z' '3d' other
% 输出：
%   tempT,tempX:    积分得到的全轨道的数据
%   T,X:            用来画图的数据，点少
%
% last modified by PH at 2013-10-12:1046
% last modified by PH at 2013-10-12:1709 修改输入支持tspan，取消对精度调节的支持
% last modified by PH at 2013-10-14:1522 重新加入对精度调节的支持，因为画HaloElliptic时不够精确了
% last modified by PH at 2013-10-24:2133 画出平动点的位置，并自动判断是 L1 还是 L2
% last modified by PH at 2013-10-25:1417 画出secxnd primary 的位置
% last modified by PH at 2013-10-04:1027 加入精度结构Tol
% last modified by PH at 2013-11-27:2108 重新编写了注释，加入 PlotType
% last modified by PH at 2014-01-06:1330 调整了输出的名字；默认画图每圈只画300个点；修改注释
% last modified by PH at 2014-10-21:1420 加入 Tol.MaxPlotPoint 的支持，控制画图的点的个数

%% 输入检测和默认参数
%%% 绘图精度
if nargin == 5 && isstruct(PlotType) % 没有输入PlotType，采用默认
    Tol = PlotType;
    PlotType = 'other';
elseif nargin == 5 && ischar(PlotType) % 没有转入精度，采用默认
    Tol.RelTol = 1e-11;
    Tol.AbsTol = 1e-11;
elseif nargin == 4 % 没有转入精度和绘图方式，采用默认
    Tol.RelTol = 1e-11;
    Tol.AbsTol = 1e-11;
    PlotType = 'other';
elseif nargin == 3
    Tol.RelTol = 1e-11;
    Tol.AbsTol = 1e-11;
    PlotType = 'other';
elseif nargin == 2
    e = 0;
    Tol.RelTol = 1e-11;
    Tol.AbsTol = 1e-11;
    PlotType = 'other';
end

%% 可调参数
if nargin < 7 && isfield(Tol,'MaxPlotPoint')
    MaxPlotPoint = Tol.MaxPlotPoint;
else
    MaxPlotPoint = 1000 * max( round( diff(tspan([1,end]))/2/pi ), 1 ); % 画图点的个数
end


%% 主程序
% ode设置
odeOptions = odeset('RelTol',Tol.RelTol,'AbsTol',Tol.AbsTol);
if isfield(Tol,'MaxStep')
    odeOptions = odeset(odeOptions, 'MaxStep', Tol.MaxStep);
end
if nargin <= 3 % 未输入tspan时，自动生成
    odeOptions = odeset(odeOptions, 'Events',@(f,X)HaloEventSecxndCross(f,X,X0));
    tspan = [0,10];
%     PlotType = 'other';
    PlotType = '3d';
end
% 积分得到轨道
[TAll, XAll] = ode45(@(t,X)HaloOde(t,X,mu,e), tspan, X0, odeOptions);
% 插值得到画图数据
TPlot = linspace(TAll(1),TAll(end),MaxPlotPoint); % 插值，只画 MaxPlotPoint 个点的图像
XPlot = interp1(TAll, XAll, TPlot); % 如果需要更高精度，可以输出 tempT 和 tempX 手动画图

%% 根据 PlotType 绘图
if nargout == 0
    switch PlotType
        case {'x-y'} % 只画 x-y
            plot(XPlot(:,1),XPlot(:,2));
            grid on; hold on; axis equal;
            plot(XPlot(1,1),XPlot(1,2),'b*');
            plot(XPlot(end,1),XPlot(end,2),'ro');
            title('x-y');
            if X0(1) < 1-mu
                plot(LibrationPoint(mu,'L1'),0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            else
                plot(LibrationPoint(mu,'L2'),0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            end
            plot(1-mu,0, 'bo','MarkerFacecolor','b','MarkerSize',7);
        case {'y-z'} % 只画 y-z
            plot(XPlot(:,2),XPlot(:,3));
            grid on; hold on; axis equal;
            plot(XPlot(1,2),XPlot(1,3),'b*');
            plot(XPlot(end,2),XPlot(end,3),'ro');
            title('y-z');
        case {'x-z'} % 只画 x-z
            plot(XPlot(:,1),XPlot(:,3),'k');
            grid on; hold on; axis equal;
%             plot(X(1,1),X(1,3),'b*');
            plot(XPlot(end,1),XPlot(end,3),'ro');
            title('x-z');
            if X0(1) < 1-mu
                plot(LibrationPoint(mu,'L1'),0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            else
                plot(LibrationPoint(mu,'L2'),0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            end
            plot(1-mu,0, 'bo','MarkerFacecolor','b','MarkerSize',7);
        case {'3d'} % 只画 3d
            plot3(XPlot(:,1),XPlot(:,2),XPlot(:,3));
            grid on; hold on; axis equal;
            plot3(XPlot(1,1),XPlot(1,2),XPlot(1,3),'b*');
            plot3(XPlot(end,1),XPlot(end,2),XPlot(end,3),'ro');
%             xlabel('x'); ylabel('y'); zlabel('z');
            if X0(1) < 1-mu
                plot3(LibrationPoint(mu,'L1'),0,0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            else
                plot3(LibrationPoint(mu,'L2'),0,0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            end
            plot3(1-mu,0,0, 'bo','MarkerFacecolor','b','MarkerSize',7);
        otherwise % 默认全画
            subplot(2,2,1);
            plot3(XPlot(:,1),XPlot(:,2),XPlot(:,3));
            grid on; hold on; axis equal;
            plot3(XPlot(1,1),XPlot(1,2),XPlot(1,3),'b*');
            plot3(XPlot(end,1),XPlot(end,2),XPlot(end,3),'ro');
            xlabel('x'); ylabel('y'); zlabel('z');
            if X0(1) < 1-mu
                plot3(LibrationPoint(mu,'L1'),0,0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            else
                plot3(LibrationPoint(mu,'L2'),0,0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            end
            plot3(1-mu,0,0, 'bo','MarkerFacecolor','b','MarkerSize',7);
            
            subplot(2,2,2);
            plot(XPlot(:,1),XPlot(:,2));
            grid on; hold on; axis equal;
            plot(XPlot(1,1),XPlot(1,2),'b*');
            plot(XPlot(end,1),XPlot(end,2),'ro');
            title('x-y');
            if X0(1) < 1-mu
                plot(LibrationPoint(mu,'L1'),0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            else
                plot(LibrationPoint(mu,'L2'),0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            end
            plot(1-mu,0, 'bo','MarkerFacecolor','b','MarkerSize',7);
            
            subplot(2,2,3);
            plot(XPlot(:,2),XPlot(:,3));
            grid on; hold on; axis equal;
            plot(XPlot(1,2),XPlot(1,3),'b*');
            plot(XPlot(end,2),XPlot(end,3),'ro');
            title('y-z');
            
            subplot(2,2,4);
            plot(XPlot(:,1),XPlot(:,3));
            grid on; hold on; axis equal;
            plot(XPlot(1,1),XPlot(1,3),'b*');
            plot(XPlot(end,1),XPlot(end,3),'ro');
            title('x-z');
            if X0(1) < 1-mu
                plot(LibrationPoint(mu,'L1'),0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            else
                plot(LibrationPoint(mu,'L2'),0, 'kx','MarkerFacecolor','c','MarkerSize',7);
            end
            plot(1-mu,0, 'bo','MarkerFacecolor','b','MarkerSize',7);
    end

end



end