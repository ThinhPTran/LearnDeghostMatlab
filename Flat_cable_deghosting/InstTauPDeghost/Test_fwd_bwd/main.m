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
fx = 0;


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



%% Input data parameters
nx3 = 3;
nt3 = 500;
ntau3 = nt3;
np3 = 1001;

dt3 = 0.004;
dx3 = 8.0;
dtau3 = dt3;
dp3 = 2.0/(vwater*(np3-1));
fp3 = -1.0/vwater;
fx3 = 0;


nf3 = nt3;
nkx3 = nx3;
df3 = 1.0/((nt3)*dt3);
dkx3 = 2.0*pi/((nx3)*dx3);


t3=(0:dt3:(nt3-1)*dt3)';
f3=(0:df3:(nf3-1)*df3)';
kx3=dkx3*[0:floor(nkx3/2) -floor(nkx3/2):-1]';
x3=(fx3:dx3:fx3+(nx3-1)*dx3)';
tau3=(0:dtau3:(ntau3-1)*dtau3)';
p3 = (fp3:dp3:-fp3)';
omega3 = 2.0*pi*f3;



%% Extend on the boundaries
npad3 = 0;
exnx3 = nx3 + 2*npad3;
exnt3 = 2*nt3;
exntau3 = exnt3;
exnp3 = np3;

exdt3 = dt3;
exdx3 = dx3;
exdtau3 = dtau3;
exdp3 = dp3;
exfp3 = fp3;
exfx3 = fx3-npad3*exdx3;

exnf3 = exnt3;
exnkx3 = exnx3;
exdf3 = 1/((exnt3)*exdt3);
exdkx3 = 2.0*pi/((exnx3)*dx3);

ext3=(0:exdt3:(exnt3-1)*exdt3)';
exf3=(0:exdf3:(exnf3-1)*exdf3)';
exkx3=exdkx3*[0:floor(exnkx3/2) -floor(exnkx3/2):-1]';
exx3=(exfx3:exdx3:exfx3+(exnx3-1)*exdx3)';
extau3=(0:exdtau3:(exntau3-1)*exdtau3)';
exxp3 = (exfp3:exdp3:-exfp3)';
exomega3 = 2.0*pi*exf3;




%% Deghost parameters
d_src = 16.0;
d_rec = 32.0;

fmin=1
fmax=120
nfb=15



%% input data
fid = fopen('data_with_ghost','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


figure();
imagesc(input);
title('input');


%% Take three traces
threetrcs=input(:,200:202);


figure();
imagesc(threetrcs);
title('three\_traces');


%% Test TauP transform

% TauP transform for the whole gather
disp('TauP transform for the whole gather');
tic
input_taup=taup_fwd(nf,f,nx,x,np,p,input);
toc


% TauP transform for three traces
disp('TauP transform for 3 traces');
tic
threetrcs_taup=taup_fwd(nf3,f3,nx3,x3,np3,p3,threetrcs);
toc


figure();
imagesc(input_taup);
title('input\_taup');


figure();
imagesc(threetrcs_taup);
title('threetrcs\_taup');



%% Test Inverse TauP transform

% Inverse TauP transform for the whole gather
disp('Inverse TauP transform for the whole gather');
tic
input_inv=taup_bwd(iter_cg,nt,t,nf,f,nx,x,np,p,input_taup);
toc


% Inverse TauP transform for three traces
disp('Inverse TauP transform for 3 traces');
tic
threetrcs_inv=taup_bwd(iter_cg,nt3,t3,nf3,f3,nx3,x3,np3,p3,threetrcs_taup);
toc


figure();
imagesc(input_inv);
title('input\_inv');


figure();
imagesc(threetrcs_inv);
title('threetrcs\_inv');









