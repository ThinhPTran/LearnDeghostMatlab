clear all
close all
clc


vwater = 1500;
amplitude = 100;


nt = 101;
dt = 0.004;
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


nfb=nt;
fmin=0.001;
res=5.0;


z=15;

tau = 2.0*z/vwater;
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;
tau1 = (index-1-floor(nt/2))*dt;



%% Get primary
noise = 0.0005*amplitude*randn(nt,1);

wlet1=rickerwave(40,0.2,t);
wlet2=rickerwave(80,0.2,t);  


primary_sig1=zeros(nt,1)+circshift(wlet1,[0 0]);
primary_sig2=zeros(nt,1)+circshift(wlet2,[0 0]); 



figure;
plot(t,primary_sig1,t,primary_sig2);
legend('primary\_sig1','primary\_sig2'); 


%% Synthetize data with ghost
withghost1 = convolution2(primary_sig1,F) ;%+ noise;




% [fspecpsig1,tmp] = analyse_spec_fwd(primary_sig1,1,nt,dt,nfb,fmin,fmax,res);
% [fspecpsig2,tmp] = analyse_spec_fwd(primary_sig2,1,nt,dt,nfb,fmin,fmax,res);
% 
% 
% 
% figure()
% imagesc(abs(fspecpsig1)); 
% 
% 
% figure();
% imagesc(abs(fspecpsig2)); 



fpsig1=fft(primary_sig1); 
fpsig2=fft(primary_sig2);
fwgh1=fft(withghost1); 



%fpsig11=abs(fwgh1).*exp(1i*angle(fpsig1));
fpsig11=smooth(abs(fwgh1)).*exp(1i*angle(fpsig1));
fpsig21=abs(fpsig1).*exp(1i*angle(fpsig2)); 


primary_sig11=real(ifft(fpsig11)); 
primary_sig21=real(ifft(fpsig21)); 


figure();
plot(t,primary_sig1,t,primary_sig11,t,withghost1);
legend('primary\_sig1','primary\_sig11','withghost1'); 


figure(); 
plot(t,primary_sig1,t,primary_sig21); 
legend('primary_sig1','primary_sig21'); 





