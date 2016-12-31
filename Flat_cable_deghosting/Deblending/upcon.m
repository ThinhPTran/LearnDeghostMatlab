clear all
close all
clc


num_z_step = 8;
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


fid = fopen('receiver_data_p','rb');
primary = fread(fid,[nt nx],'single');
fclose(fid);

primary_fk = fft2(primary);

fid = fopen('primary_fk.bin','wb');
fwrite(fid,abs(primary_fk),'single');
fclose(fid);

fid = fopen('primary_fk_angle.bin','wb');
fwrite(fid,angle(primary_fk),'single');
fclose(fid);


wavefield = upward_extr(primary, nf, omega, nkx, kx, num_z_step, dz);
wavefield_fk = fft2(wavefield);


fid = fopen('upcon.bin','wb');
fwrite(fid,wavefield,'single');
fclose(fid);

fid = fopen('upcon_fk.bin','wb');
fwrite(fid,abs(wavefield_fk),'single');
fclose(fid);

fid = fopen('upcon_fk_angle.bin','wb');
fwrite(fid,angle(wavefield_fk),'single');
fclose(fid);







