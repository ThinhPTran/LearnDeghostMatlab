close all
clear all
clc

iter_cg = 10000;
epsilon = 0.5;

nt = 2001;
n2 = 1001;
dt = 0.004;
tmax = dt*(nt-1);
df = 1./tmax;
fmax = df*(nt-1);

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;

fmin=0;
fmax=100;

fminindex = floor(fmin/df) + 1;
fmaxindex = floor(fmax/df) + 1;

z = 50.0;
tau = 2.0*z/1500;
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;

fid = fopen('primary.bin','r');
synthetic = fread(fid,[nt n2],'float');
fclose(fid);


%% Get primary
primary_sig = synthetic(1:nt,floor(n2/2));
noise = randn(2001,1)/10.0;
primary_sig = primary_sig + 1./cosh(100*(t-1)) +  1./cosh(100*(t-2)) + 1./cosh(100*(t-4)) + 1./cosh(100*(t-5)) + 1./cosh(100*(t-6));
figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F);% + noise;
figure(2)
plot(t,withghost);


%% Frequency domain
fwithghost = fft(withghost);
figure(3)
plot(f,abs(fwithghost));


%% Phase 
figure(4)
anglefwithghost = angle(fwithghost);
smoothanglefwithghost = smooth(f,anglefwithghost,0.05,'loess');
plot(f,angle(fwithghost),'yellow');
hold on
plot(f,smoothanglefwithghost,'red');
hold off




%% De-noise 
% fwithghost = fft(withghost);
% 
% for i=fmaxindex+1:nt
%   fwithghost(i)=0;
% end
% 
% for i = nt:-1:floor(nt/2)
%   fwithghost(i) = conj(fwithghost(nt - i + 2));
% end
% 
% withghost = real(ifft(fwithghost));
% 
% tmp = forward_PEF(withghost,iter_cg,epsilon);
% forwardfilter = zeros(nt,1);
% forwardfilter(floor(nt/2) + 1) = tmp(1);
% forwardfilter(floor(nt/2) + 2) = tmp(2);
% 
% forwardpredictwithghost = convolution2(withghost,forwardfilter);
% forwardpredictwithghost = circshift(forwardpredictwithghost,[1 0]);
% 
% tmp = backward_PEF(withghost,iter_cg,epsilon);
% backwardfilter = zeros(nt,1);
% backwardfilter(floor(nt/2) + 1) = tmp(1);
% backwardfilter(floor(nt/2) + 2) = tmp(2);
% 
% backwardpredictwithghost = convolution2(withghost,backwardfilter);
% backwardpredictwithghost = circshift(backwardpredictwithghost,[-2 0]);
% 
% withghost = (forwardpredictwithghost + backwardpredictwithghost)/2.0;
% 
% 
% figure(3);
% plot(t,withghost);
% 
% 

%% Create Mirror data
mirror_data = -withghost;
figure(5)
plot(t,mirror_data);




%% Find the Primary (Deghosting) 

tic;

   
 F1 = cal_filter(withghost, -withghost, iter_cg, epsilon);


figure(6)
plot(t,F1);




%% Create ghost filter

maxind = 1;
maxvalue = F1(1);
for i = 2:nt
    if ( F1(i) > maxvalue ) 
       maxindex = i;
       maxvalue = F1(i);
    end
end

if (maxindex < (floor(nt/2) + 1)) 
  tmp = floor(nt/2) + 1 - maxindex;
  maxindex = floor(nt/2) + 1 + tmp;
end

F1 = zeros(nt,1);
F1(floor(nt/2)+1) = 1.0;
F1(maxindex) = -1.0;

figure(7);
plot(t,F1);


%% Calculate P

% P = cal_primary(F1, withghost, iter_cg, epsilon);
%  
%  
%  figure(7);
%  plot(t,P);
 
 
 %% Calculate P second way.
 tau2 = (maxindex - floor(nt/2) - 1)*dt
 %tau2 = 0.064;
 minustshift = exp(1i*tau2*omega);
 a = 1.0;
 F2 = 1 - a*minustshift;
 denominator = 1 + a*a -  2*a*cos(omega*tau2) + 1e-10;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P2 = real(ifft(fP));
 figure(8)
 plot(t,P2);
 

toc;