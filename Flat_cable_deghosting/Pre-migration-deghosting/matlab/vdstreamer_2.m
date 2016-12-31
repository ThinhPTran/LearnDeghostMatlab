close all
clear all
clc


nt = 1001;
dt = 0.004;
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


z = 24;
tau = 2.0*z/1500;
amplitude = 100;
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = amplitude*randn(nt,1)*0.05;
wlet = ricker(40,0.004);
primary_sig =  zeros(nt,1) +amplitude*1./cosh(100*(t-0.3));
for i = 1:13
   primary_sig(i+10) = 1.0*amplitude*wlet(i); 
   primary_sig(i+150) = -0.8*amplitude*wlet(i);
   primary_sig(i+300) = 1.3*amplitude*wlet(i);
   primary_sig(i+450) = 0.5*amplitude*wlet(i);
   primary_sig(i+600) = 1.1*amplitude*wlet(i);
   primary_sig(i+750) = 0.3*amplitude*wlet(i);
   primary_sig(i+900) = -1.3*amplitude*wlet(i);
end

figure(1)
plot(t,primary_sig);


addback = 0.1;
ntwindow = nt;
dtwindow = dt;
dfwindow = 1./((ntwindow-1)*dtwindow);
tmaxwindow = (ntwindow-1)*dtwindow;
fmaxwindow = (ntwindow-1)*dfwindow;

twindow = (0:dtwindow:tmaxwindow)';
fwindow = (0:dfwindow:fmaxwindow)';
omegawindow = 2*pi*fwindow;


zmax = 75;
zmin = 15;
tmaxcheck = 1.5*zmax*2/1500;
fmaxcheck = fmaxwindow/4;
fmincheck = fix(1/tmaxcheck);
tmincheck1 = 1./fmaxcheck;
tmincheck2 = 2*zmin/1500;

tmincheck = max([tmincheck1 tmincheck2]);
tmaxcheck;
fmincheck;
fmaxcheck;

ntcheck = fix(tmaxcheck/dtwindow) + 1;
%ntcheck = nt;
tcheck = 0:dtwindow:(ntcheck-1)*dtwindow;
testmean = zeros(ntcheck,1);
maxvalue = zeros(ntcheck,1);
check = zeros(ntcheck,1);


denominatorstore = zeros(ntwindow,ntcheck);
F2store = zeros(ntwindow,ntcheck);
denominatore = zeros(ntwindow,1);
F2 = zeros(ntwindow,1);
withghost = zeros(ntwindow,1);
fwithghost = zeros(ntwindow,1);
absfwithghost = zeros(ntwindow,1);
fP = zeros(ntwindow,1);
P = zeros(ntwindow,1);
G = zeros(ntwindow,1);


%% Synthetize data with ghost
tmp = convolution2(primary_sig,F);%  + noise;
withghost = tmp;
figure(2)
plot(t,withghost);


tic


for i = 1:ntcheck
    tmptau = (i-1)*dtwindow;
    denominatorstore(:,i) = 2 -  2*cos(omegawindow*tmptau) + 2.0;
    minustshift = exp(1i*tmptau*omegawindow);
    F2store(:,i) = 1 - minustshift;
end


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);
figure(3)
plot(f,absfwithghost);


average = mean(abs(withghost));
maxP = max(withghost);


for iter = 1:ntcheck
 F2 = F2store(:,iter);
 denominator = denominatorstore(:,iter);
 fP =( fwithghost.*F2)./denominator;
 for i = ntwindow:-1:floor(ntwindow/2)
  fP(i) = conj(fP(ntwindow - i + 2));
 end

 P = real(ifft(fP));
 G = withghost - P;
 
 maxvalue(iter) = max(P);
 test = G + circshift(P,[iter-1 0]);
 testmean(iter) = 1./mean(abs(test)/maxP);

end


for i = 1:ntcheck
    check(i) = testmean(i)+maxvalue(i);
end


[tmp ind] = max(check);


%% Recover
 tau2 = (ind -1)*dtwindow;
 F2 = F2store(:,ind);
 denominator = 2 -  2*cos(omegawindow*tau2) + 0.2;
 fP =( fwithghost.*F2)./denominator;
 for i = ntwindow:-1:floor(ntwindow/2)
   fP(i) = conj(fP(ntwindow - i + 2));
 end

 P = real(ifft(fP));
 G = withghost-P;

 
 P = (addback*withghost + P)/(1.0+addback);
 fP = fft(P);
 absfP = abs(fP);
 
 
 toc
 
 
 tau
 tau2
 
 
 figure(4)
 plot(twindow,withghost,'blue',twindow,P,'red',twindow,G,'green');
 
 
 figure(5);
 plot(tcheck,testmean,'red',tcheck,maxvalue,'green',tcheck,check,'blue');
 
 
 figure(6);
 plot(twindow,absfwithghost,'red',twindow,absfP,'blue');




