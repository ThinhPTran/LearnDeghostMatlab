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


z = 6;
fpeak = 20;
tau = 2.0*z/1500;
amplitude = 100;
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = amplitude*randn(nt,1)*0.05;
wlet = ricker(fpeak,dt);
primary_sig =  zeros(nt,1) +amplitude*1./cosh(100*(t-0.3));
for i = 1:13
  primary_sig(i+10) = 1.0*amplitude*wlet(i); 
   primary_sig(i+150) = -1.0*amplitude*wlet(i);
   primary_sig(i+300) = 1.0*amplitude*wlet(i);
   primary_sig(i+450) = 1.0*amplitude*wlet(i);
   primary_sig(i+600) = 1.0*amplitude*wlet(i);
   primary_sig(i+750) = 1.0*amplitude*wlet(i);
   primary_sig(i+900) = -1.0*amplitude*wlet(i);
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
tmp = convolution2(primary_sig,F)  + noise;
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
 
 maxvalue(iter) = max(P)/average;
 test = G + circshift(P,[iter-1 0]);
 testmean(iter) = 1./mean(abs(test)/average);

end


for i = 1:ntcheck
    check(i) = 1.5*testmean(i)+maxvalue(i);
end


[tmp ind] = max(check);

%% Recover
 tau2 = (ind -1)*dtwindow; 
 
 
 %% Calculate P, 1st way
 F2 = F2store(:,ind);
 denominator = 2 -  2*cos(omegawindow*tau2) + 0.00001;
 fP =( fwithghost.*F2)./denominator;
 for i = ntwindow:-1:floor(ntwindow/2)
   fP(i) = conj(fP(ntwindow - i + 2));
 end

 
 %% Calculate P, 2nd way
 threshold = 1.0
 Fshift = 1 - exp(-1i*tau2*omegawindow);
 absFshift = abs(Fshift);
 
 for i = 1:ntwindow
    if (absFshift(i) < threshold) 
      Fshift(i) = sqrt(0.5*threshold*threshold)+1i*sqrt(0.5*threshold*threshold);
    end
 end
 
 absFshift = abs(Fshift);
 fP1 = fwithghost./Fshift;
  for i = ntwindow:-1:floor(ntwindow/2)
   fP1(i) = conj(fP1(ntwindow - i + 2));
 end
 
 
 P = real(ifft(fP));
 P1= real(ifft(fP1));
 G = withghost-P;
 G1 = withghost-P1;

 
 %% Calculate P 
 weight = abs(P+G);
 maxweight = max(weight);
 normweight = weight/maxweight;
 P2 = P.*normweight;
 
 
 P = (addback*withghost + P)/(1.0+addback);
 P1 = (addback*withghost + P1)/(1.0+addback);
 P2 = (addback*withghost + P2)/(1.0+addback);
 
 
 fP = fft(P);
 fP1 = fft(P1);
 fP2 = fft(P2);
 absfP = abs(fP);
 absfP1 = abs(fP1);
 absfP2 = abs(fP2);
 
 
 max1 = max(absfwithghost);
 max2 = max(absfP2);
 scale = max1/max2;
 
 fP = fP*scale;
 fP1 = fP1*scale;
 fP2 = fP2*scale;
 absfP = abs(fP);
 absfP1 = abs(fP1);
 absfP2 = abs(fP2);
 
 P = real(ifft(fP));
 P1 = real(ifft(fP1));
 P2 = real(ifft(fP2));
 
 toc
 
 
 tau
 tau2
 
 
 figure(4)
 plot(twindow,withghost,'blue',twindow,P2,'red');
 
 
 figure(5);
 plot(tcheck,testmean,'red',tcheck,maxvalue,'green',tcheck,check,'blue');
 
 
 figure(6);
 plot(twindow,absfwithghost,'red',twindow,absfP2,'blue');
 
 
