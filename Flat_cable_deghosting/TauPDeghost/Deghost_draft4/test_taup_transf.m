clear all
close all
clc



%% Control parameters
vwater = 1500;
recover_eps = 0.05;
p_control = 1.0;
ctr_n = 1.0;
iter_cg = 1000;
res = 5.0;



%% Input data parameters
nx = 401;
nt = 500;
ntau = nt;
np = 1001;

dt = 0.004;
dx = 8.0;
dtau = dt;
dp = 2.0/(vwater*(np-1));
fp = -1.0/vwater;
fx = -1600;


nf = nt;
nkx = nx;
df = 1.0/((nt)*dt);
dkx = 2.0*pi/((nx)*dx);


t=(0:dt:(nt-1)*dt)';
f=(0:df:(nf-1)*df)';
kx=dkx*[0:floor(nkx/2) -floor(nkx/2):-1]';
x=(fx:dx:fx+(nx-1)*dx)';
tau=(0:dtau:(ntau-1)*dtau)';
p = (fp:dp:-fp)';
omega = 2.0*pi*f;



%% Extend on the boundaries
npad = 100;
exnx = nx + 2*npad;
exnt = 2*nt;
exntau = exnt;
exnp = np;

exdt = dt;
exdx = dx;
exdtau = dtau;
exdp = dp;
exfp = fp;
exfx = fx-npad*exdx;

exnf = exnt;
exnkx = exnx;
exdf = 1/((exnt)*exdt);
exdkx = 2.0*pi/((exnx)*dx);

ext=(0:exdt:(exnt-1)*exdt)';
exf=(0:exdf:(exnf-1)*exdf)';
exkx=exdkx*[0:floor(exnkx/2) -floor(exnkx/2):-1]';
exx=(exfx:exdx:exfx+(exnx-1)*exdx)';
extau=(0:exdtau:(exntau-1)*exdtau)';
exxp = (exfp:exdp:-exfp)';
exomega = 2.0*pi*exf;



%% Deghost parameters
d_src = 16.0;
d_rec = 32.0;

fmin=1
fmax=120
nfb=4



%% input data
fid = fopen('data_with_ghost','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


%% TauP transform
disp('TauP transform');
tic
input_taup = taup_fwd(nf,f,nx,x,np,p,input);
toc

fid = fopen('fwd_taup.bin','wb');
fwrite(fid,input_taup,'single');
fclose(fid);


%% Inverse TauP transform
disp('Inverse TauP transform');
tic
inv_taup = taup_bwd(iter_cg,nt,t,nf,f,nx,x,np,p,input_taup);
toc

fid = fopen('inv_taup.bin','wb');
fwrite(fid,inv_taup,'single');
fclose(fid);


sum_in=sum(input,2);
sum_inv=sum(inv_taup,2);


fsum_in=abs(fft(sum_in));
fsum_inv=abs(fft(sum_inv));


figure();
plot(f,fsum_in,f,fsum_inv);
legend('fsum\_in','fsum\_inv');


dbfsum_in=20*log(fsum_in)-140;
dbfsum_inv=20*log(fsum_inv)-140;


figure();
plot(f,dbfsum_in,f,dbfsum_inv);
legend('dbfsum\_in','dbfsum\_inv');









