function rho = RTBPPrimaryDistance(e0,f)
% 计算 RTBP 中主天体间的距离
%
%   created by PH at 2017-04-12:1859
%   last modified by PH at 2017-09-18:1002 更正错误：不应该使用输入的 f0，坐标定义为近地点时 f mod 2pi = 0

rho = (1-e0.^2) ./ (1+e0.*cos(f));

end
    