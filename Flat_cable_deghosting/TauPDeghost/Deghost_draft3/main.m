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
nfb=25



%% input data
fid = fopen('data_with_ghost','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


input_fk = fft2(input);

fid = fopen('data_with_ghost_fk','wb');
fwrite(fid,abs(input_fk),'single');
fclose(fid);


exinput = zeros(exnt, exnx);
exinput(1:nt,npad+1:npad+nx) = input;
 


%% TauP transform
tic
exinput_taup = taup_fwd(exnf,exf,exnx,exx,exnp,exxp,exinput);
toc

fid = fopen('receiver_data_with_ghost_n1_500_taup.bin','wb');
fwrite(fid,exinput_taup,'single');
fclose(fid);



%% Deghost in TauP domain
tic
exinput_taup_deghost=taup_deghost(p_control,vwater,d_src,d_rec,recover_eps,exnt,exdt,ext,exomega,exnp,exdp,exxp,exinput_taup);
toc


exerr_taup = abs(exinput_taup_deghost - exinput_taup);

fid = fopen('exerr_taup.bin','wb');
fwrite(fid,exerr_taup,'single');
fclose(fid);


%% Clean up
tic
exinput_taup_deghost_clean=cleanup(exinput_taup,exinput_taup_deghost,exnp,exnt,exdt,nfb,fmin,fmax,res);
toc


exerr_taup_clean = abs(exinput_taup_deghost_clean - exinput_taup);

fid = fopen('exerr_taup_clean.bin','wb');
fwrite(fid,exerr_taup_clean,'single');
fclose(fid);


fid = fopen('receiver_data_with_ghost_n1_500_taup_deghost.bin','wb');
fwrite(fid,exinput_taup_deghost,'single');
fclose(fid);

fid = fopen('receiver_data_with_ghost_n1_500_taup_deghost_clean.bin','wb');
fwrite(fid,exinput_taup_deghost_clean,'single');
fclose(fid);



%% Inverse TauP transform
tic
extmp_inv = taup_bwd(iter_cg,exnt,ext,exnf,exf,exnx,exx,exnp,exxp,exinput_taup_deghost_clean);
toc

tmp_exoutput = clean_non_physical_fk(ctr_n,vwater,exnt,ext,exnf,exf,exomega,exnx,exx,exnp,exxp,exnkx,exkx,extmp_inv);


output = tmp_exoutput(1:nt,npad+1:npad+nx);
% output = extmp_inv(1:nt,npad+1:npad+nx);
output_fk = fft2(output);



fid = fopen('receiver_data_with_ghost_n1_500_invtaup_deghost.bin','wb');
fwrite(fid,output,'single');
fclose(fid);


fid = fopen('receiver_data_with_ghost_n1_500_invtaup_deghost_fk.bin','wb');
fwrite(fid,abs(output_fk),'single');
fclose(fid);


mean2(abs(output))









