clear all
close all
clc



%% Control parameters
vwater=1500;
recover_eps=0.5;
p_control=1.0;
ctr_n=1.1;
iter_cg=1000;
res=5.0;



%% Input data parameters
nx=101;
nt=500;
ntau=nt;
np=1001;

dt=0.004;
dx=8.0;
dtau=dt;
dp=2.0/(vwater*(np+1));
fp=-1.0/vwater+dp;
fx=100;


nf=nt;
nkx=nx;
df=1.0/((nt)*dt);
dkx=2.0*pi/((nx)*dx);


t=(0:dt:(nt-1)*dt)';
f=(0:df:(nf-1)*df)';
kx=dkx*[0:floor(nkx/2) -floor(nkx/2):-1]';
x=(fx:dx:fx+(nx-1)*dx)';
tau=(0:dtau:(ntau-1)*dtau)';
p=(fp:dp:-fp)';
omega=2.0*pi*f;



%% Extend on the boundaries
npad=50;
exnx=nx+2*npad;
exnt=2*nt;
exntau=exnt;
exnp=np;

exdt=dt;
exdx=dx;
exdtau=dtau;
exdp=dp;
exfp=fp;
exfx=fx-npad*exdx;

exnf=exnt;
exnkx=exnx;
exdf=1/((exnt)*exdt);
exdkx=2.0*pi/((exnx)*dx);

ext=(0:exdt:(exnt-1)*exdt)';
exf=(0:exdf:(exnf-1)*exdf)';
exkx=exdkx*[0:floor(exnkx/2) -floor(exnkx/2):-1]';
exx=(exfx:exdx:exfx+(exnx-1)*exdx)';
extau=(0:exdtau:(exntau-1)*exdtau)';
exxp=(exfp:exdp:-exfp)';
exomega=2.0*pi*exf;



%% Deghost parameters
x0=-20;
x_noff=100;
alpha=atan(1.0/8.0);

d=(x+x0)*tan(alpha);

fmin=1;
fmax=120;
nfb=4;

delta_d=60;



%% input data
fid=fopen('slant_with_ghost_case1.bin','r');
input=fread(fid,[nt nx],'single');
fclose(fid);


fid=fopen('flat_realcable_case1.bin','r');
flatinput=fread(fid,[nt, nx],'single');
fclose(fid);



input_taup=taup_fwd(nf,f,nx,x,np,p,input);
flatinput_taup=taup_fwd(nf,f,nx,x,np,p,flatinput);


tic
input_taup_dg=taup_slant_dg(recover_eps,alpha,vwater,nt,dt,np,dp,input_taup);
toc



figure();
imagesc(input_taup);
title('input\_taup');


figure();
imagesc(flatinput_taup);
title('flatinput\_taup');


figure();
imagesc(input_taup_dg);
title('input\_taup\_dg');



input_inv=taup_bwd(iter_cg,nt,t,nf,f,nx,x,np,p,input_taup);


%% Output one t trace
ttrace=input_taup(100,:);
fid=fopen('ttrace.bin','wb');
fwrite(fid,ttrace,'single');
fclose(fid);





















