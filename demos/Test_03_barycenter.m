clear;

%% convert to earth

cspice_kclear;
miceRootFolder = '/Users/GroupMacBai/Codes/Tools/SPICE/generic_kernels/';
cspice_furnsh([miceRootFolder 'lsk/naif0011.tls']);
cspice_furnsh([miceRootFolder 'spk/planets/de432s.bsp']);
cspice_furnsh([miceRootFolder 'pck/gm_de431.tpc']);

% nondimensionanl --> dimensional
% after conversion, assumed in EMBR mice frame
%%% epoches
et0 = cspice_str2et('2017/11/13 00:00:00'); % [s] % the ephemeris epoch corresponding to initialEpoches(1)
et1 = cspice_str2et('2037/11/13 00:00:00');
et = linspace(et0,et1,100001); % [s] % the time interval matters
%%% states
earth = cspice_spkezr('Earth',et,'J2000','none','Earth');
moon = cspice_spkezr('Moon',et,'J2000','none','Earth');
bc = cspice_spkezr('Earth Moon Barycenter',et,'J2000','none','Earth');
%%%
GMearth = cspice_bodvrd( 'earth', 'GM', 1 ); % [km^3/s^2]
GMmoon = cspice_bodvrd( 'moon', 'GM', 1 ); % [km^3/s^2]
mu = GMmoon / (GMearth + GMmoon);
bc2 = (1-mu)*earth + mu*moon;

% difference
disp(max(max(abs(bc-bc2))));

%%
% figure(93);
% clf;
% plot(bc(1,:),'ro'); hold on;
% plot(bc2(1,:),'b*'); hold on;