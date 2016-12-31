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
z = 15;
fpeak = 40;
zmin_check = z*0.3;
zmax_check = z*1.5;

tmin_check = 2.0*zmin_check/vwater;
tmax_check = 2.0*zmax_check/vwater;

dt_check = 0.001;
t_check = (tmin_check:dt_check:tmax_check)';
nt_check = length(t_check);


tau = 2.0*z/vwater
index = floor(tau/dt) + floor(nt/2) + 1;


F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;
tau1 = (index-1-floor(nt/2))*dt;



%% Get primary
noise = 0.0001*amplitude*randn(nt,1);
% wlet = ricker(fpeak,dt);
wlet = rickerwave(fpeak,0.2,t);
primary_sig = zeros(nt,1);% +amplitude*1./cosh(100*(t-0.3));


primary_sig = primary_sig + circshift(wlet,[0 0]);
% primary_sig = primary_sig + circshift(wlet,[50 0]);


figure;
plot(t,primary_sig);


%% Synthetize data with ghost
withghost = convolution2(primary_sig,F) + noise;


figure
plot(t,withghost);


normtestmean = zeros(nt_check,1);
normmaxvalue = zeros(nt_check,1);
normderiv = zeros(nt_check,1);
normderiv2 = zeros(nt_check,1);
normcheck = zeros(nt_check,1);
normtest = zeros(nt_check,1);


average = mean(abs(withghost));


for iter = 15:15
    
 newP = zeros(size(withghost));
 P = zeros(size(withghost));

 tau2 = t_check(iter);
 
 for j_iter = 1:30
    
   P = newP; 
   fP = fft(P);

   minustshift = exp(1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2) + recover_eps;
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
 iter
 tau2
 fP = fft(P);
 shiftfP = fP.*exp(-1i*tau2*omega);
 shiftP = real(ifft(shiftfP));
 test = shiftP+G;
 normtestmean(iter) = average/mean(abs(test));
 
 % Calculate normcheck
 normcheck(iter) = 1.0/mean(abs(P));
 
 
 figure();
 plot(t,P,t,shiftP,t,G,t,withghost,t,test);
 legend('P','shiftP','G','withghost','test');
 
 
end


%% Normalize
normmaxvalue=normmaxvalue/(max(abs(normmaxvalue)));
normtestmean=normtestmean/(max(abs(normtestmean)));
normcheck=normcheck/(max(abs(normcheck)));

% Calculate normderiv
normderiv=abs(circshift(normmaxvalue,[1 0]) - circshift(normmaxvalue,[-1 0]));
normderiv2=abs(circshift(normmaxvalue,[1 0]) - normmaxvalue + circshift(normmaxvalue,[-1 0]));
for iter=1:nt_check
  if (normderiv(iter)>0) 
    normderiv(iter)=1.0/normderiv(iter); 
  end
end
% normderiv=1.0./normderiv;
normderiv=normderiv/(max(abs(normderiv)));
normderiv2=normderiv2/(max(abs(normderiv2)));


normtest=(normmaxvalue+normderiv2+normtestmean)/3.0;


figure()
plot(t_check,normderiv2,t_check,normmaxvalue,t_check,normtestmean,t_check,normtest);
legend('normderiv2','normmaxvalue','normtestmean','normtest');
% plot(t_check,normcheck,t_check,normmaxvalue,t_check,normtestmean,t_check,normderiv,t_check,normderiv2,t_check,normtest);
% legend('normcheck','normmaxvalue','normtestmean','normderiv','normderiv2','normtest');


[tmp ind] = max(normtest);


%% Recover
 tau
 tau1
 tau2 = (ind-1)*dt_check+tmin_check
 tau2 = tau1
 
 
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
 
 end
 
figure
plot(t,withghost,'blue',t,P,'red',t,primary_sig,'green');
legend('withghost','Primary recovered','Primary');


















