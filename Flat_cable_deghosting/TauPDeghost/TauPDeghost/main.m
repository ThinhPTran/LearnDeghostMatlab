clear all
close all
clc


%% Control parameters
vwater = 1500;
recover_eps = 0.01;
p_control = 0.8;
iter_cg = 1000;



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



%% input data
fid = fopen('taup_imag_after_deghost','r');
input_taup_deghost = fread(fid,[exnt np],'single');
fclose(fid);



%% Inverse TauP transform
tic
tmp_inv = taup_bwd(iter_cg,exnt,ext,exnf,exf,nx,x,np,p,input_taup_deghost);
toc

tmp_output = tmp_inv; %taup_clean(1.0,vwater,exnt,ext,exnf,exf,exomega,nx,x,np,p,nkx,kx,tmp_inv);


output = tmp_output(1:nt,:);


fid = fopen('receiver_data_with_ghost_n1_500_invtaup_deghost.bin','wb');
fwrite(fid,output,'single');
fclose(fid);


mean2(abs(output))









