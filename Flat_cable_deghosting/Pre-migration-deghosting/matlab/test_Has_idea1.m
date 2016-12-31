close all
clear all
clc


iter_cg = 10000;
finding_eps = 2.0;
recover_eps = 2.0;
epsilon = 5.0;
vwater = 1500;
amplitude = 100;


nt = 101;
n2 = 101;
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


z = 15;
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
noise = 0.0005*amplitude*randn(nt,1);
% wlet = ricker(fpeak,dt);
wlet = rickerwave(fpeak,0.2,t);
primary_sig = zeros(nt,1);% +amplitude*1./cosh(100*(t-0.3));


primary_sig = primary_sig + circshift(wlet,[0 0]);
% primary_sig = primary_sig + circshift(wlet,[50 0]);


figure;
plot(t,primary_sig);


%% Synthetize data with ghost
withghost = convolution2(primary_sig,F) ;%+ noise;


figure
plot(t,withghost);


testmean = zeros(nt_check,1);
maxvalue = zeros(nt_check,1);
check = zeros(nt_check,1);


average = mean(abs(withghost));


for iter = ntcheck_min:ntcheck_max
    
 newG = zeros(size(withghost));
 G = zeros(size(withghost));

 tau2 = (iter-1)*dt_check;
 
 for j_iter = 1:30
    
   G = newG; 
   fG = fft(G);

   minustshift = exp(-1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2) + recover_eps;
   fwithghost = fft(withghost);
   fG =( fwithghost.*F2 + recover_eps*fG)./denominator;
   for j = nt:-1:floor(nt/2)
    fG(j) = conj(fG(nt - j + 2));
   end
 

   newG = real(ifft(fG));
   diff = abs(newG - G);
   
   f = mean(abs(diff))/mean(abs(newG));
   
   if (f < 0.01) 
      break; 
   end
 
 end
 
 j_iter;
 
 P = withghost - G;

 
 maxvalue(iter) = max(G)/average;
 test = G + circshift(P,[iter-1 0]);
 testmean(iter) = average/mean(abs(test));
 check(iter) = 1.0/mean(abs(P));
 
 
end


figure()
plot(t_check,check);


[tmp ind] = max(check);


%% Recover
 tau
 tau1
 tau2 = (ind-1)*dt_check
 tau2 = 0.02
 
 
 P = zeros(size(withghost));
 
 for j_iter = 1:30
    
   fP = fft(P);

   minustshift = exp(1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2) + recover_eps;
   fwithghost = fft(withghost);
   fP =( fwithghost.*F2 + recover_eps*fP)./denominator;
   for j = nt:-1:floor(nt/2)
    fP(j) = conj(fP(nt - j + 2));
   end
 

   P = real(ifft(fP));
   
%    figure()
%    plot(t,P);
 
 end
 
 
 G = zeros(size(withghost));
 
 for j_iter = 1:30
    
   fG = fft(G);

   minustshift = exp(-1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2) + recover_eps;
   fwithghost = fft(withghost);
   fG =( fwithghost.*F2 + recover_eps*fG)./denominator;
   for j = nt:-1:floor(nt/2)
    fG(j) = conj(fG(nt - j + 2));
   end
 

   G = real(ifft(fG));

   
   
 end

 
 
figure
plot(t,withghost,'blue',t,P,'red',t,primary_sig,'green');
legend('withghost','Primary recovered','Primary');

figure
plot(t,withghost,'blue',t,G,'red',t,primary_sig,'green');
legend('withghost','Ghost recovered','Primary');