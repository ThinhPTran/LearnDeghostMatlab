close all
clear all
clc


iter_cg = 10000;
finding_eps = 2.0;
recover_eps = 0.0001;
epsilon = 0.05;
vwater = 1500;
amplitude = 100;


nt = 1001;
n2 = 1001;
dt = 0.004;
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;



%% t for checking
dt_check = 0.001;
t_check = (0:dt_check:tmax)';
nt_check = length(t_check);


z = 4;
fpeak = 40;
zmin_check = 3;
zmax_check = 80;


tau = 2.0*z/vwater;
taumin_check = 2.0*zmin_check/vwater;
taumax_check = 2.0*zmax_check/vwater;
index = floor(tau/dt) + floor(nt/2) + 1;

ntcheck_min = floor(taumin_check/dt_check)+1;
ntcheck_max = floor(taumax_check/dt_check)+1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;
tau1 = (index-1-floor(nt/2))*dt;



%% Get primary
noise = 0.05*amplitude*randn(nt,1);
% wlet = ricker(fpeak,dt);
wlet = rickerwave(fpeak,0.5,t);
primary_sig = zeros(nt,1);% +amplitude*1./cosh(100*(t-0.3));


primary_sig = primary_sig + circshift(wlet,[10-1 0]);


figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F);%  + noise;


figure(2)
plot(t,withghost);


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);
figure(3)
plot(f,absfwithghost);


check = zeros(nt_check,1);


for iter = ntcheck_min:ntcheck_max

 tau2 = (iter-1)*dt_check;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
%  denominator = ones(size(minustshift));
 denominator = 2 -  2*cos(omega*tau2) + finding_eps;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P = real(ifft(fP));
 G = withghost - P;
 
 test = G + circshift(P,[iter-1 0]);
 
 average = mean(abs(P));
 check(iter) = max(P)/average;
 
end


figure(4);
plot(t_check,check,'blue');
legend('check');


[tmp ind] = max(check);



%% Recover
 tau
 tau1
 tau2 = (ind-1)*dt_check
%  tau2 = tau
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2) + recover_eps;
%  denominator = ones(size(minustshift));
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end
 

 P2 = real(ifft(fP));
 G = withghost - P2;
 
 
figure(5)
plot(t,withghost,'blue',t,P2,'red',t,primary_sig,'green');
legend('withghost','Primary recovered','Primary');




