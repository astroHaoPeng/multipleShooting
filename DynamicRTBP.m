function Xdot = RTBPOde(f,X,mu,e,number_vectors,delta_f)
% Restricted Three-Body Problem 动力学方程
%
%   以真近点角为独立变量；CRTBP 时即为时间，ERTBP 时为角度
%
%   调用格式：
%       Xdot = RTBPOde(f,X,mu,e,number_vectors,delta_f)
%
%   Input:
%       f,X,mu,e:   
%       number_vectors: 向量化计算时 X 中包含的点数
%       delta_f:        向量化的各段与输入的时刻 f 之间的时间差：向量化要求各段的时间相等，f+delta_f 得到各段的真实时间
%
%   Output:
%       Xdot:       与输入的 X 相同维数的向量
%
%
%   注意：勿轻易修改，小的改动会导致小的数值误差，先前计算的轨道库的精度就会受到影响
%
% 	created by PH at 2014-04-21:1842: 从 HaloOde 修改而来，以后都调用这个
% 	last modified by PH at 2014-08-12:1247 完善注释；增加对向量输入的支持（未验证）
%   PH at 2014-11-06:1436 需要对多输入进行验证

%% 输入检测
if nargin < 4
    error('Not enougth input! Provide input at least as (f,X,mu,e).');
elseif nargin < 5
	if length(X) > 6 % ~exist('number_vectors','var')
		error('Please input number of vectorized states: number_vectors as the fifth input.');
	elseif length(X) == 6
		number_vectors = 1;
		delta_f = 0;
	end
elseif nargin == 6
	if length(X)/6 ~= number_vectors
		error('Size of X and value of number_vectors do not fit.');
	end
	if length(delta_f)/6 ~= number_vectors
		error('Size of delta_f and value of number_vectors do not fit.');
	end
end

%% 向量化时恢复输入为矩阵
% 状态按列存储
f = f + reshape(delta_f,[],1);
X_mat = reshape(X,6,number_vectors);

%% 由 symbolic.m 推导来的ODE
x = X_mat(1,:);
y = X_mat(2,:);
z = X_mat(3,:);
dx = X_mat(4,:);
dy = X_mat(5,:);
dz = X_mat(6,:);

r2cubic = ((mu + x - 1).^2 + y.^2 + z.^2).^(3/2);
r1cubic = ((mu + x).^2 + y.^2 + z.^2).^(3/2);
kappa = (e.*cos(f) + 1);

Xdot = [
    dx;
    dy;
    dz;
     2*dy + (x + ((2*mu + 2*x).*(mu - 1))./(2*r1cubic) - (mu.*(2*mu + 2*x - 2))./(2*r2cubic))./kappa;
    -2*dx + (y - (mu*y)./r2cubic + (y.*(mu - 1))./r1cubic)./kappa;
    -z + (z - (mu*z)./r2cubic + (z.*(mu - 1))./r1cubic)./kappa];

Xdot = reshape(Xdot, size(X));

end