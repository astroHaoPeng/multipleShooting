syms GM
syms x y z vx vy vz
syms px py pz
syms phix phiy phiz phivx phivy phivz

X = [x,y,z].'; % S/C
P = [px,py,pz].'; % planet


v = [vx, vy, vz].';
a = - GM / sum( (X-P).^2 )^(3/2) * (X-P) - GM / sum( P.^2 )^(3/2) * P

dX = [v; a];

phiX = [phix,phiy,phiz,phivx,phivy,phivz].';

simplify( jacobian(dX, [x,y,z,vx,vy,vz]) * phiX )