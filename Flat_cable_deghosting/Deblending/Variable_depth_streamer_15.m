close all
clear all
clc

iter_cg = 10000;
epsilon = 0.1;
vwater = 1500;
amplitude = 1;


nt = 1001;
n2 = 1001;
dt = 0.004;
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


z = 80;
zmin_check = 5;
zmax_check = 200;


tau = 2.0*z/vwater;
taumin_check = 2.0*zmin_check/vwater;
taumax_check = 2.0*zmax_check/vwater;
index = floor(tau/dt) + floor(nt/2) + 1;

ntcheck_min = floor(taumin_check/dt)+1;
ntcheck_max = floor(taumax_check/dt)+1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = 0.05*amplitude*randn(nt,1);
wlet = ricker(40,0.004);
primary_sig =  zeros(nt,1) +amplitude*1./cosh(100*(t-0.3));
for i = 1:13
   primary_sig(i+10) = 1.0*amplitude*wlet(i); 
   primary_sig(i+150) = 0.8*amplitude*wlet(i);
   primary_sig(i+300) = 1.3*amplitude*wlet(i);
   primary_sig(i+450) = 0.5*amplitude*wlet(i);
   primary_sig(i+600) = 1.1*amplitude*wlet(i);
   primary_sig(i+750) = 0.3*amplitude*wlet(i);
   primary_sig(i+900) = 1.0*amplitude*wlet(i);
end

figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F)  + noise;
%withghost = noise;
figure(2)
plot(t,withghost);


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);
figure(3)
plot(f,absfwithghost);


testmean = zeros(nt,1);
maxvalue = zeros(nt,1);
check = zeros(nt,1);


average = mean(abs(withghost));
averagespec = mean(absfwithghost);


for iter = ntcheck_min:ntcheck_max

 tau2 = (iter-1)*dt;
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
 
 maxvalue(iter) = max(P2)/average;
 test = G2 + circshift(P2,[iter-1 0]);
 testmean(iter) = 1./mean(abs(test)/average);

end


for i = 1:nt
    check(i) = testmean(i)+maxvalue(i);
end



figure(4);
plot(t,check,'blue',t,maxvalue,'black',t,testmean,'red');
legend('check','maxvalue','testmean');


[tmp ind] = max(check);

%% Recover
 tau 
 tau2 = (ind -1)*dt
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2);
 stable_p = 0.1*mean(abs(denominator));
 denominator = denominator + 0.1;
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