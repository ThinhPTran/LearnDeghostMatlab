close all
clear all
clc



recover_eps = 2.0;
vwater = 1500;
amplitude = 1.0;
s2n = 10000;


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
z = 10;
fpeak = 80;
zmin_check = z*0.7;
zmax_check = z*1.3;

dt_check = 0.0001;
tmin_check = floor((2.0*zmin_check/vwater)/dt_check)*dt_check;
tmax_check = floor((2.0*zmax_check/vwater)/dt_check)*dt_check;


t_check = (tmin_check:dt_check:tmax_check)';
nt_check = length(t_check);


tau = 2.0*z/vwater
index = floor(tau/dt) + floor(nt/2) + 1;


F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;
%tau1 = (index-1-floor(nt/2))*dt



%% Get primary
noise = (1.0/s2n)*amplitude*randn(nt,1);
% wlet = ricker(fpeak,dt);
wlet = rickerwave(fpeak,0.2,t);
primary_sig = zeros(nt,1) +amplitude*1./cosh(800*(t-0.2));


%primary_sig = zeros(nt,1) + circshift(wlet,[0 0]);
%primary_sig = primary_sig + circshift(wlet,[10 0]);



figure;
plot(t,primary_sig);
legend('primary_sig');


% figure;
% plot(f,fprimary_sig);
% legend('fprimary_sig');


%% Synthetize data with ghost
%withghost = convolution2(primary_sig,F) + noise;
fprimary_sig = fft(primary_sig);
fwithghost = fprimary_sig.*(1-exp(-1i*tau*omega));
for j = nt:-1:floor(nt/2)-1
    fwithghost(j) = conj(fwithghost(nt - j + 2));
end

% figure();
% plot(f,abs(fwithghost));
% legend('fwithghost');

withghost = real(ifft(fwithghost))+noise+(1.0/s2n)*amplitude*1./cosh(800*(t-0.3));


figure
plot(t,withghost);
legend('withghost');

tic
%% Calculate norms
normtestmean = zeros(nt_check,1);
normmaxvalue = zeros(nt_check,1);
normcheck = zeros(nt_check,1);
normtest = zeros(nt_check,1);


average = mean(abs(withghost));


for iter = 1:nt_check
    
 newP = zeros(size(withghost));
 P = zeros(size(withghost));

 ctau = t_check(iter);
 
 for j_iter = 1:30
    
   P = newP; 
   fP = fft(P);

   minustshift = exp(1i*ctau*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*ctau) + recover_eps;
   fwithghost = fft(withghost);
   fP =( fwithghost.*F2 + recover_eps*fP)./denominator;
   for j = nt:-1:floor(nt/2)
    fP(j) = conj(fP(nt - j + 2));
   end
 

   newP = real(ifft(fP));
   diff = abs(newP - P);
   
   f = mean(abs(diff))/mean(abs(newP));
   
   if (f < 0.01) 
      break; 
   end
 
 end
 
 j_iter;
 
 G = withghost - P;
 

 % Calculate normmax value
 normmaxvalue(iter) = max(P)/average;
 
 % Calculate normtestmean value
 fP = fft(P);
 shiftfP = fP.*exp(-1i*ctau*omega);
 for j = nt:-1:floor(nt/2)
    shiftfP(j) = conj(shiftfP(nt - j + 2));
 end
 shiftP = real(ifft(shiftfP));
 test = shiftP+G;
 normtestmean(iter) = average/mean(abs(test));
 
 % Calculate normcheck
 normcheck(iter) = 1.0/mean(abs(P));
 
 
%  figure();
%  plot(t,P,t,shiftP,t,G,t,withghost,t,test);
%  legend('P','shiftP','G','withghost','test');
 
 
end


% Normalize
normmaxvalue=normmaxvalue/(max(abs(normmaxvalue)));
normtestmean=normtestmean/(max(abs(normtestmean)));
normcheck=normcheck/(max(abs(normcheck)));



% Calculate norm for test
% if (z<=6)
% normtest=(normmaxvalue+2.0*normderiv2+2.0*normtestmean+normcheck)/7.0;
% else 
normtest=(normmaxvalue+normtestmean+normcheck)/3.0;    
% end


figure()
plot(t_check,normmaxvalue,t_check,normtestmean,t_check,normtest,t_check,normcheck);
legend('normmaxvalue','normtestmean','normtest','normcheck');


[tmp ind] = max(normtest);



%% Recover
 est_tau = (ind-1)*dt_check+tmin_check
 
 
 P = zeros(size(withghost));
 
 for j_iter = 1:30
    
   fP = fft(P);

   minustshift = exp(1i*est_tau*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*est_tau) + recover_eps;
   fwithghost = fft(withghost);
   fP =(fwithghost.*F2 + recover_eps*fP)./denominator;
   for j = nt:-1:floor(nt/2)
    fP(j) = conj(fP(nt - j + 2));
   end
 

   P = real(ifft(fP));
 
 end
 
toc
 
figure
plot(t,withghost,'blue',t,P,'red',t,primary_sig,'green');
legend('withghost','Primary recovered','Primary');


















