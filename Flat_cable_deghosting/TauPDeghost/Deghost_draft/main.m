clear all
close all
clc


%% Control parameters
vwater = 1500;
recover_eps = 0.01;
p_control = 1.0;
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



%% Extend parameters
exnt = 2*nt;
exnf = exnt;
exntau = exnt;
nkx = nx;
exdt = dt;
exdtau = dtau;
exdf = 1.0/((exnt)*exdt);
dkx = 2.0*pi/((nx)*dx);


ext=(0:exdt:(exnt-1)*exdt)';
exf=(0:exdf:(exnf-1)*exdf)';
kx=dkx*[0:floor(nkx/2) -floor(nkx/2):-1]';
x=(fx:dx:fx+(nx-1)*dx)';
extau=(0:exdtau:(exntau-1)*exdtau)';
p = (fp:dp:-fp)';
exomega = 2.0*pi*exf;



%% Deghost parameters
d_src = 16.0;
d_rec = 32.0;

fmin=1
fmax=120
nfb=400



%% input data
fid = fopen('data_with_ghost','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


input_fk = fft2(input);

fid = fopen('data_with_ghost_fk','wb');
fwrite(fid,abs(input_fk),'single');
fclose(fid);


exinput = [input; zeros(nt,nx)];
 


%% TauP transform
tic
exinput_taup = taup_fwd(exnf,exf,nx,x,np,p,exinput);
toc

fid = fopen('receiver_data_with_ghost_n1_500_taup.bin','wb');
fwrite(fid,exinput_taup,'single');
fclose(fid);



%% Deghost in TauP domain
tic
exinput_taup_deghost = taup_deghost(p_control,vwater,d_src,d_rec,recover_eps,exnt,exdt,ext,exomega,np,dp,p,exinput_taup);
toc



%% Clean up
exinput_taup_deghost_clean=cleanup(exinput_taup,exinput_taup_deghost,np,exnt,exdt,nfb,fmin,fmax,res);


fid = fopen('receiver_data_with_ghost_n1_500_taup_deghost.bin','wb');
fwrite(fid,exinput_taup_deghost,'single');
fclose(fid);

fid = fopen('receiver_data_with_ghost_n1_500_taup_deghost_clean.bin','wb');
fwrite(fid,exinput_taup_deghost_clean,'single');
fclose(fid);



%% Inverse TauP transform
tic
extmp_inv = taup_bwd(iter_cg,exnt,ext,exnf,exf,nx,x,np,p,exinput_taup_deghost_clean);
toc

tmp_exoutput = clean_non_physical_fk(1.15,vwater,exnt,ext,exnf,exf,exomega,nx,x,np,p,nkx,kx,extmp_inv,exinput);


output = tmp_exoutput(1:nt,:);
output_fk = fft2(output);



fid = fopen('receiver_data_with_ghost_n1_500_invtaup_deghost.bin','wb');
fwrite(fid,output,'single');
fclose(fid);


fid = fopen('receiver_data_with_ghost_n1_500_invtaup_deghost_fk.bin','wb');
fwrite(fid,abs(output_fk),'single');
fclose(fid);


mean2(abs(output))









