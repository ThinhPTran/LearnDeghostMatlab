close all
clear all
clc


iter_cg = 10000;
epsilon = 0.1;
vwater = 1500;
amplitude = 1;


nt = 500;
n2 = 1001;
dt = 0.004;
df = 1./((nt)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;



%% t for checking
dt_check = 0.001;
t_check = (0:dt_check:tmax)';
nt_check = length(t_check);


zmin_check = 5;
zmax_check = 100;


taumin_check = 2.0*zmin_check/vwater;
taumax_check = 2.0*zmax_check/vwater;


ntcheck_min = floor(taumin_check/dt_check)+1;
ntcheck_max = floor(taumax_check/dt_check)+1;




%% Get primary
fid = fopen('receiver_data_with_ghost_n1_500_taup.bin','rb');
data = fread(fid,[nt n2],'single');
fclose(fid);


withghost = data(:,500);


figure(1)
plot(t,withghost);


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);
figure(2)
plot(f,absfwithghost);


check = zeros(nt_check,1);


for iter = ntcheck_min:ntcheck_max

 tau2 = (iter-1)*dt_check;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2);
 denominator = denominator + 2.0;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P2 = real(ifft(fP));
 G2 = withghost - P2;
 
 test = G2 + circshift(P2,[iter-1 0]);
 
 check(iter) = max(P2)./mean(abs(test));
 
end


figure(3);
plot(t_check,check,'blue');
legend('check');


[tmp ind] = max(check);



%% Recover
 tau_src = 0.0213
 tau_rec = 0.0426
 tau2 = (ind -1)*dt_check
 %tau2 = 0.0426
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2);
 denominator = denominator + 0.1;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end
 

 P2 = real(ifft(fP));
 G = withghost - P2;
 
figure(4)
plot(t,withghost,'blue',t,P2,'red');
legend('withghost','Primary recovered');




