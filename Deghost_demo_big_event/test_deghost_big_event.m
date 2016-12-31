clc
clear all;
close all;


nt=2251;
nx=636;

dt = 4e-3;
hnf=floor(nt/2);
t=[0:dt:(nt-1)*dt];
df=1./(nt*dt);
f=[0:hnf -hnf:-1]*df;
omega=2.0*pi*f';

eps=0.5


z1=6.0;
z2=8.0; 
vwater=1500.0; 
tau1=2.0*z1/vwater;
tau2=2.0*z2/vwater; 



fid=fopen('datain.bin','r'); 
datain=fread(fid,[2251 636],'single'); 
fclose(fid); 


onetrc=datain(:,160); 


figure(); 
plot(t(1500:1560),onetrc(1500:1560)); 
title('onetrc'); 


finput=fft(onetrc); 


nominator1=1-exp(1i*omega*tau1);
denominator1=2.0-2.0*cos(omega*tau1)+eps;

nominator2=1-exp(1i*omega*tau2); 
denominator2=2.0-2.0*cos(omega*tau2)+eps;


foutput1=finput.*nominator1./denominator1;
foutput2=foutput1.*nominator2./denominator2; 


output=real(ifft(foutput2)); 


figure(); 
plot(t(1500:1560), onetrc(1500:1560), t(1500:1560), output(1500:1560))
legend('input','deghost'); 



figure(); 
plot(t, onetrc, t, output); 













