function dYOrigin = RTBPOdeSTM(f,Y,mu,e)
%% ERTBP 和 CRTBP 下状态转移矩阵 Phi 满足的 ODE 方程
% dY = HaloTransitionOde(f,Y,mu,e)
%   Y(1:6) 为 Phi 的状态
%   Y(7:12) 为轨道的状态
%
% 大量简化了自动推导出的公式的计算量，考虑手动推导，进一步减少计算量
%
% created by PH at 2013-12-28:1114 从之前的各种计算 Phi 的函数中独立出来的程序，方便简化修改
% last modified by PH at 2014-01-06:1316 整理注释
% last modified by PH at 2014-01-09:1556 加入了对多列 Phi 的直接支持，不需要再使用外部循环
% last modified by PH at 2014-02-22:1322 进一步精简了 dY 计算，减小了大量乘法计算，已验证无误
%
%   created by PH at 2014-09-22:1411 从 HaloTransitionOde 修改而来
%   last modified by PH at 2017-04-12:2032 comments only: 脑子抽抽啦……为什么当时把轨道状态放最后啊……


N = length(Y)/6 - 1; % Y 中包含的 Phi 的列数
YOrigin = Y;
dYOrigin = zeros(size(YOrigin));

for ii = 1:N
    Y = [ YOrigin(1+(ii-1)*6:ii*6); YOrigin(end-5:end) ]; % 从中恢复为一对一对的形式，方便循环
    
    xx = Y(1);
    yy = Y(2);
    zz = Y(3);
    dxx = Y(4);
    dyy = Y(5);
    dzz = Y(6);
    x = Y(7);
    y = Y(8);
    z = Y(9);
    % dx = Y(10);
    % dy = Y(11);
    % dz = Y(12);
    
    r1cubic = ((mu + x)^2 + y^2 + z^2)^(3/2);
    r2cubic = ((mu + x - 1)^2 + y^2 + z^2)^(3/2);
    r1power5 = ((mu + x)^2 + y^2 + z^2)^(5/2);
    r2power5 = ((mu + x - 1)^2 + y^2 + z^2)^(5/2);
    kappa = (e*cos(f) + 1);
    
    % 进一步精简乘法计算前的正确公式
%     dY = [dxx
%         dyy
%         dzz
%         2*dyy + (xx*((mu - 1)/r1cubic - mu/r2cubic + (3*mu*(2*mu + 2*x - 2)^2)/(4*r2power5) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*r1power5) + 1))/kappa + (yy*((3*mu*y*(2*mu + 2*x - 2))/(2*r2power5) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*r1power5)))/kappa + (zz*((3*mu*z*(2*mu + 2*x - 2))/(2*r2power5) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*r1power5)))/kappa
%         (yy*((mu - 1)/r1cubic - mu/r2cubic - (3*y^2*(mu - 1))/r1power5 + (3*mu*y^2)/r2power5 + 1))/kappa - 2*dxx - (zz*((3*y*z*(mu - 1))/r1power5 - (3*mu*y*z)/r2power5))/kappa + (xx*((3*mu*y*(2*mu + 2*x - 2))/(2*r2power5) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*r1power5)))/kappa
%         zz*(((mu - 1)/r1cubic - mu/r2cubic - (3*z^2*(mu - 1))/r1power5 + (3*mu*z^2)/r2power5 + 1)/kappa - 1) - (yy*((3*y*z*(mu - 1))/r1power5 - (3*mu*y*z)/r2power5))/kappa + (xx*((3*mu*z*(2*mu + 2*x - 2))/(2*r2power5) - (3*z*(2*mu + 2*x)*(mu - 1))/(2*r1power5)))/kappa
%         reshape(HaloOde(f,Y(7:12),mu,e), 6, 1)];
    
    muminus1 = mu - 1;
    two_x = 2 * x;
    two_mu = 2 * mu;
    three_mu = 3 * mu;
    
    dY = [dxx
        dyy
        dzz
        2*dyy + (xx*(muminus1/r1cubic - mu/r2cubic + (three_mu*(two_mu + two_x - 2)^2)/(4*r2power5) - (3*(two_mu + two_x)^2*muminus1)/(4*r1power5) + 1))/kappa + (yy*((three_mu*y*(two_mu + two_x - 2))/(2*r2power5) - (3*y*(two_mu + two_x)*muminus1)/(2*r1power5)))/kappa + (zz*((three_mu*z*(two_mu + two_x - 2))/(2*r2power5) - (3*z*(two_mu + two_x)*muminus1)/(2*r1power5)))/kappa
        (yy*(muminus1/r1cubic - mu/r2cubic - (3*y^2*muminus1)/r1power5 + (three_mu*y^2)/r2power5 + 1))/kappa - 2*dxx - (zz*((3*y*z*muminus1)/r1power5 - (three_mu*y*z)/r2power5))/kappa + (xx*((three_mu*y*(two_mu + two_x - 2))/(2*r2power5) - (3*y*(two_mu + two_x)*muminus1)/(2*r1power5)))/kappa
        zz*((muminus1/r1cubic - mu/r2cubic - (3*z^2*muminus1)/r1power5 + (three_mu*z^2)/r2power5 + 1)/kappa - 1) - (yy*((3*y*z*muminus1)/r1power5 - (three_mu*y*z)/r2power5))/kappa + (xx*((three_mu*z*(two_mu + two_x - 2))/(2*r2power5) - (3*z*(two_mu + two_x)*muminus1)/(2*r1power5)))/kappa
        reshape(DynamicRTBP(f,Y(7:12),mu,e), 6, 1)];
%     if norm(dY-dYY) ~= 0
%         disp(norm(dY-dYY)); % 验证是否相等
%     end
    
    
    dYOrigin(1+(ii-1)*6:ii*6) = reshape(dY(1:6),6,1);
    
end

dYOrigin(end-5:end) = dY(end-5:end);

end