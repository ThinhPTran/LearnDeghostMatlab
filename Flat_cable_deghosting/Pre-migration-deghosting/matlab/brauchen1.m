close all
clear all
clc


iter_cg = 10000;
finding_eps = 2.0;
recover_eps = 0.2;
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
withghost = convolution2(primary_sig,F) + noise;


figure
plot(t,withghost);


checkerror1 = zeros(nt_check,1);
checkerror2 = zeros(nt_check,1);
checkerror3 = zeros(nt_check,1);
checkerror4 = zeros(nt_check,1);
check = zeros(nt_check,1);


average = mean(abs(withghost));


%for iter = ntcheck_min:ntcheck_max
for iter = 20:20

 newP = zeros(size(withghost));
 P = zeros(size(withghost));
 newG = zeros(size(withghost));
 G = zeros(size(withghost));


 tau2 = (iter-1)*dt_check;
 
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
 
 P = newP;
 
 
 for j_iter = 1:30
    
   G = newG; 
   fG = fft(G);

   minustshift = exp(-1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2) + recover_eps;
   fwithghost = fft(withghost);
   fG =( fwithghost.*F2 + recover_eps*fG)./denominator;
   for j = nt:-1:floor(nt/2)
    fG(j) = conj(fP(nt - j + 2));
   end
 

   newG = real(ifft(fG));
   diff = abs(newG - G);
   
   f = mean(abs(diff))/mean(abs(newG));
   
   if (f < 0.01) 
      break; 
   end
 
 end
 
 G = newG;
 
 
 fP = fft(P);
 minustshift = exp(-1i*tau2*omega);
 fG1 = -fP.*minustshift;
 G1 = real(ifft(fG1));

 
 fG = fft(G);
 minustshift = exp(1i*tau2*omega);
 fP1 = -fG.*minustshift;
 P1 = real(ifft(fP1));
 
 
 checkerror1(iter) = mean(abs(withghost - (P+G)));
 checkerror2(iter) = mean(abs(withghost - (P1+G1)));
 checkerror3(iter) = mean(abs(G1-G));
 checkerror4(iter) = mean(abs(P1-P));
 
end


checkerror = checkerror1+checkerror2+checkerror3+checkerror4;

figure()
plot(t_check,checkerror1,t_check,checkerror2);%,t_check,checkerror3,t_check,checkerror4);
legend('checkerror1','checkerror2');%,'checkerror3','checkerror4');



derror = circshift(checkerror,1)-circshift(checkerror,-1);


figure()
plot(t_check,1./abs(derror));
legend('derror');


[tmp ind] = max(check);


figure()
plot(t,P,t,G);
legend('P','G');

figure()
plot(t,P1,t,G1);
legend('P1','G1');


%% Recover
%  tau
%  tau1
%  tau2 = (ind-1)*dt_check
%  
%  
%  P = zeros(size(withghost));
%  
%  for j_iter = 1:30
%     
%    fP = fft(P);
% 
%    minustshift = exp(1i*tau2*omega);
%    F2 = 1 - minustshift;
%    denominator = 2 -  2*cos(omega*tau2) + recover_eps;
%    fwithghost = fft(withghost);
%    fP =( fwithghost.*F2 + recover_eps*fP)./denominator;
%    for j = nt:-1:floor(nt/2)
%     fP(j) = conj(fP(nt - j + 2));
%    end
%  
% 
%    P = real(ifft(fP));
%  
%  end
%  
% figure
% plot(t,withghost,'blue',t,P,'red',t,primary_sig,'green');
% legend('withghost','Primary recovered','Primary');