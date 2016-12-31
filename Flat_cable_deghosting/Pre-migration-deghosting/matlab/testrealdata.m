close all
clear all
clc

iter_cg = 10000;
epsilon = 0.001;

nt = 4000
n2 = 576
dt = 0.002
tmax = dt*(nt-1)
dx = 25
xmax = dx*(n2-1)
df = 1./(tmax)
dkx = 1./xmax
fmax = (nt-1)*df
kxmax = (n2-1)*dkx

t = (0:dt:tmax)';
f = (0:df:fmax)';
x = (0:dx:xmax)';
kx = (0:dkx:kxmax)';
omega = 2*pi*f;


%% Get data
data = zeros(nt,n2);

fid = fopen('input.bin','r');
data = fread(fid,[nt n2],'float');
fclose(fid);


fxdata = fft2(data);

[tmpff tmpkk] = meshgrid(f,kx);
ff = tmpff';
kk = tmpkk';
mesh(ff,kk,abs(fxdata));
colormap hsv
colorbar


