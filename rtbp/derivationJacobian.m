syms mu
syms x y z vx vy vz
syms px py pz
syms phix phiy phiz phivx phivy phivz

X = [x,y,z].'; % S/C
V = [vx, vy, vz].';

r1 = sqrt( (x+mu)^2 + y^2 + z^2 );
r2 = sqrt( (x-1+mu)^2 + y^2 + z^2 );

U = 1/2*(x^2+y^2) + (1-mu)/r1 + mu/r2 + mu*(1-mu);


A = -[diff(U,x); diff(U,y); diff(U,z)]

dX = [V; A];

phiX = [phix,phiy,phiz,phivx,phivy,phivz].';

simplify( jacobian(dX, [x,y,z,vx,vy,vz]) * phiX )