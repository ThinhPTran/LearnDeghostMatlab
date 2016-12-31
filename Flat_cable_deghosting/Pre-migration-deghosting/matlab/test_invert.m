close all
clear all
clc

iter_cg = 10000;
epsilon = 0.5;


dt = 0.004;
T = 8;
df = 1/T;
nt = floor(T/dt) + 1
n2 = 1001;

t = (0:dt:T)';
f = (0:df:(nt-1)*df)';
omega = 2*pi*f;

fmin=0;
fmax=100;

fminindex = floor(fmin/df) + 1;
fmaxindex = floor(fmax/df) + 1;

z = 50.0;
tau = 2.0*z/1500
tau = 0.064
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;



%% Get primary
noise = randn(2001,1)/10.0;
primary_sig = 1./cosh(100*(t-2)) + 1./cosh(100*(t-4));
figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F) + noise;
figure(2)
plot(t,withghost);


%% First method
% tau2 = 0.064;
% F1 = zeros(nt,1);
% F1(floor(nt/2)+1) = 1.0;
% F1(floor(nt/2)+1 + floor(tau2/dt) ) = -1.0;
% 
% figure(6);
% plot(t,F1);
% 
% 
% %% Calculate P
% 
% P = cal_primary(F1, withghost, iter_cg, epsilon);
%  
%  
%  figure(7);
%  plot(t,P);


%% Second method
 tau2 = 0.064;
 minustshift = exp(1i*tau2*omega);
 a = 1.0;
 F2 = 1 - a*minustshift;
 
 avg = mean(abs(F2))
 
 denominator = 1 + a*a -  2*a*cos(omega*tau2) + avg;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P2 = real(ifft(fP));
 figure(3)
 plot(t,P2);
 
 