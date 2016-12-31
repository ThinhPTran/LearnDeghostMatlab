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
wlet2=rickerwave(60,0.2,t);  


primary_sig1=zeros(nt,1)+circshift(wlet1,[0 0]);
primary_sig2=zeros(nt,1)+circshift(wlet2,[0 0]); 



figure;
plot(t,primary_sig1,t,primary_sig2);
legend('primary\_sig1','primary\_sig2'); 


fsig1=fft(primary_sig1); 
fsig2=fft(primary_sig2); 

absfsig1=abs(fsig1); 
absfsig2=abs(fsig2); 

dbabsfsig1=mag2db(absfsig1); 
dbabsfsig2=mag2db(absfsig2); 


figure(); 
plot(f,absfsig1,f,absfsig2);
legend('absfsig1','absfsig2'); 


figure(); 
plot(f,dbabsfsig1,f,dbabsfsig2); 
legend('dbabsfsig1','dbabsfsig2'); 














