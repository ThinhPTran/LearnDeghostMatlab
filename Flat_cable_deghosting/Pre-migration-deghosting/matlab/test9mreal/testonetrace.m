close all
clear all
clc

iter_cg = 10000;
epsilon = 0.1;
stable_p = 1.0;
addback = 0.1;

nt = 4000;
n2 = 576;
dt = 0.002;
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;



%% Get primary
fid = fopen('9mghost.bin','r');
data = fread(fid,[nt n2],'single');
withghost = data(1401:1500,250);


nt = 100;
n2 = 576;
dt = 0.002;
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


figure(1)
plot(t,withghost);


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);
figure(2)
plot(f,absfwithghost);


zmax = 50;
zmin = 3;
tmaxcheck = 1.0*zmax*2/1500;
fmaxcheck = fmax/4;
fmincheck = fix(1/tmaxcheck);
tmincheck1 = 1./fmaxcheck;
tmincheck2 = 2*zmin/1500;

tmincheck = 2*zmin/1500;
tmaxcheck;
fmincheck;
fmaxcheck;

ntcheck = fix(tmaxcheck/dt) + 1;
tcheck = 0:dt:(ntcheck-1)*dt;
testmean = zeros(ntcheck,1);
maxvalue = zeros(ntcheck,1);
baseonf = zeros(ntcheck,1);
check = zeros(ntcheck,1);

findcheck = fix(tmincheck/dt) + 1;
lindcheck = fix(tmaxcheck/dt) + 1;


maxP = max(withghost);
average = mean(abs(withghost));
averagespec = mean(absfwithghost);


for iter = 1:lindcheck

 tau2 = (iter-1)*dt;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2) + 2.0;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P = real(ifft(fP));
 
 maxvalue(iter) = max(P)/average;

end


figure(3);
plot(tcheck,maxvalue);


for iter = 1:lindcheck

 tau2 = (iter-1)*dt;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2) + 2.0;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P = real(ifft(fP));
 G2 = withghost - P;
 
 test = (G2 + circshift(P,[iter-1 0]))/average;
 testmean(iter) = 1./mean(abs(test));

end


figure(4);
plot(tcheck,testmean);



for i = 1:ntcheck
    check(i) = testmean(i)+maxvalue(i);
end



figure(5);
plot(tcheck,check,'blue',tcheck,maxvalue,'black',tcheck,testmean,'red');


[tmp ind] = max(check);

%% Recover
 tau2 = (ind -1)*dt
 if ( abs(tau2 - 0.012) < 0.04 )
     
   minustshift = exp(1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2) + 0.0000000000000001;
   fwithghost = fft(withghost);
   fP =( fwithghost.*F2)./denominator;
   for i = nt:-1:floor(nt/2)
    fP(i) = conj(fP(nt - i + 2));
   end
   
   P = real(ifft(fP));
 
   G = withghost - P;
   
 else 
     
   P = withghost;
   G = withghost - P;
   
 end
 
 %% Calculate P 
 weight = abs(P+G);
 maxweight = max(weight);
 normweight = weight/maxweight;
 P = P.*normweight;
 
 
 fP = fft(P);
 absfP = abs(fP);
 
 
 %% scale energy
 max1 = max(absfwithghost);
 max2 = max(absfP);
 scale = max1/max2;
 fP = fP*scale;
 absfP = abs(fP);

 
 P = real(ifft(fP));
 
 
 P = (addback*withghost + P)/(1.0+addback);
 
 
 toc
 
 
 tau2
 
 
 meancheck = mean(abs(check));
 maxcheck = max(abs(check));
 noisecriteria = maxcheck/meancheck
 
 
 figure(6)
 plot(t,withghost,'blue',t,P,'red',t,G,'green');
 
 
 figure(7);
 plot(tcheck,testmean,'red',tcheck,maxvalue,'green',tcheck,check,'blue');
 
 
 figure(8);
 plot(t,absfwithghost,'red',t,absfP,'blue');  











