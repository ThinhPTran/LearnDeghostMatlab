clear all
close all
clc


sigma=sqrt(0.00004);


dt=0.0004;
tmax=2.0;
t=0:dt:tmax;
Nt=length(t);

df=1.0/tmax;
fmax=df*(Nt-1);
f=0:df:fmax;
Nf=Nt;



y1=exp(-0.5*((t-0.04).^2)/(sigma^2));
fy1=fft(y1); 
y2=(1-((t-0.04).^2)/(sigma^2)).*exp(-0.5*((t-0.04).^2)/(sigma^2));
fy2=fft(y2);


figure(1);
plot(t,y1);
legend('Gaussian');


figure(2);
plot(f,abs(fy1));
legend('Fourier of Gaussian');


figure(3);
plot(t,y2);
legend('Ricker Wavelet');


figure(4);
plot(f,abs(fy2));
legend('Fourier of Ricker Wavelet');


% ind = find(t > 0.95 & t < 1.05);
% wavelet = y2(ind);
% 
% length(ind)
% 
% 
% figure(5);
% plot(wavelet);

fid = fopen('source_wavelet','w');
fwrite(fid,y1,'single');
fclose(fid);







