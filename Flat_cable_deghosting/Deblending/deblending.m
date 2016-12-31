% deblending
clear all
close all
clc


epsilon = 100;
iter_cg = 1000;
dz = 8.0;
z = 0.0;


nt = 500;
dt = 0.004;
tmin = 0.0;
tmax = (nt-1)*dt;
t = (tmin:dt:tmax)';

nf = nt;
domega = 2.0*pi/(nt*dt);
omega = [0:floor(nt/2) -floor(nt/2)+1:-1]'*domega;


nx = 401;
dx = 8.0;
xmin = 400.0;
xmax = (nx-1)*dx;
x = xmin:dx:xmax;

nkx = nx;
dkx = 2.0*pi/(nx*dx);
kxmin = 0.0;
kxmax = (nkx-1)*dkx;
kx = [0:floor(nx/2) -floor(nx/2):-1]*dkx;

tic

fid = fopen('data_with_ghost','rb');
data = fread(fid,[nt nx],'single');
fclose(fid);


G = -data;
P = data;


ghost_model = - upward_extr(P, nf, omega, nkx, kx, 8, dz);


fid = fopen('ghost_model','wb');
fwrite(fid,ghost_model,'single');
fclose(fid);

estimated_primary = leastsub(data, ghost_model, nx, nt, epsilon, iter_cg);

fid = fopen('estimated_primary','wb');
fwrite(fid,estimated_primary,'single');
fclose(fid);

toc















