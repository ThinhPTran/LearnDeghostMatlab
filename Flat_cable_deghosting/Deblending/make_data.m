clear all
close all
clc


%% Input data
vwater = 1500;


%% Input data parameters
nx = 401;
nt = 500;
ntau = nt;
np = 1001;

dt = 0.004;
dx = 8.0;
dtau = dt;
dp = 2.2/(vwater*(np-1));
fp = -1.1/vwater;
fx = -1600;


nf = nt;
nkx = nx;
df = 1.0/((nt)*dt);
dkx = 2.0*pi/((nx)*dx);


t=(0:dt:(nt-1)*dt)';
f=(0:df:(nf-1)*df)';
kx=dkx*[0:floor(nkx/2) -floor(nkx/2):-1]';
x=(fx:dx:-fx)';
tau=(0:dtau:(ntau-1)*dtau)';
p = (fp:dp:-fp)';
omega = 2.0*pi*f;


fid = fopen('receiver_data_p','rb');
primary = fread(fid,[nt nx],'single');
fclose(fid);


fid = fopen('receiver_data_g','rb');
ghost = fread(fid,[nt nx],'single');
fclose(fid);


%% Remove unwanted signal
% delta_t = 6.5e-4*delta_x

for i_iter = 1:nt
  for j_iter = 1:nx
    if ((t(i_iter) < 6.4e-4*x(j_iter) + 0.12) || (t(i_iter) < -6.4e-4*x(j_iter) + 0.12))
      primary(i_iter, j_iter) = 0.0;
      ghost(i_iter,j_iter) = 0.0;
    end
  end
end


fid = fopen('receiver_data_p_after_removing','wb');
fwrite(fid,primary,'single');
fclose(fid);

fid = fopen('receiver_data_g_after_removing','wb');
fwrite(fid,ghost,'single');
fclose(fid);



ghost = -ghost;


data_with_ghost = primary + ghost;


fid = fopen('data_with_ghost','wb');
fwrite(fid,data_with_ghost,'single');
fclose(fid);


data_with_ghost_fk = fft2(data_with_ghost);


fid = fopen('data_with_ghost_fk','wb');
fwrite(fid,abs(data_with_ghost_fk),'single');
fclose(fid);


primary_fk = fft2(primary);
fprimary = fft(primary);


fid = fopen('primary_fk','wb');
fwrite(fid,abs(primary_fk),'single');
fclose(fid);


fid = fopen('fprimary','wb');
fwrite(fid,abs(fprimary),'single');
fclose(fid);















