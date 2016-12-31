clear all
close all
clc



%% Control parameters
vwater=1500;
recover_eps=0.15;
p_control=1.0;
ctr_n=1.1;
iter_cg=1000;
res=5.0;



%% Input data parameters
nx=240;
nt=501;
ntau=nt;
np=2001;

dt=0.004;
dx=12.5;
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


%% Prepare a wavelet
sig=sinc(100.0*(t-1.0));


% plot to see
figure();
plot(t,sig);
title('sig');



%% Generate ghosted signal
tau=0.1;
fsig=fft(sig);
foutsig=fsig.*(1.0 - exp(-1i.*omega*tau));
outsig=real(ifft(foutsig));


% Plot to see
figure();
plot(t,outsig);
title('outsig');






























