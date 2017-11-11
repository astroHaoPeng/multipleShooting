function Test_01_propagator()
mu = 3.040229489168351e-06;
X0 = [0.994057479631466   0.000000000000000  -0.004296775120051   0.000000000000000  -0.017253209942360  -0.000000000000000];

odeOptions = odeset('AbsTol',1e-13,'RelTol',1e-13,'MaxStep',0.001) ;

tic;
[~,PVHalo45] = ode45(...
    @(f,X)DynamicRTBP(f,X,mu,0),...
    linspace(0,50*2*pi,10001),...
    X0,...
    odeOptions);
toc;

tic;
[~,PVHalo113] = ode113(...
    @(f,X)DynamicRTBP(f,X,mu,0),...
    linspace(0,50*2*pi,10001),...
    X0,...
    odeOptions);
toc;

%%
figure(92);
clf;
plot(linspace(0,100*HaloPeriod(0,X0,mu),10001)/pi,PVHalo45(:,1),'.'); hold on;
plot(linspace(0,100*HaloPeriod(0,X0,mu),10001)/pi,PVHalo113(:,1),'x'); hold on;
% plot(linspace(0,50*2*pi,10001)/2/pi,PVHalo45(:,1)-PVHalo113(:,1),'.'); hold on;

max(max(abs(PVHalo45-PVHalo113)))


end





function Xdot = DynamicRTBP(f,X,mu,e)
%% 由 symbolic.m 推导来的ODE
x = X(1,:);
y = X(2,:);
z = X(3,:);
dx = X(4,:);
dy = X(5,:);
dz = X(6,:);

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