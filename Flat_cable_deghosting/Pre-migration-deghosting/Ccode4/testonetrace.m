close all
clear all
clc


vwater = 5000;

iter_cg = 10000;
epsilon = 0.001;
stable_p = 1.0;
addback = 0.0;

nt = 3001;
n2 = 100;
dt = 0.00066666;
df = 1./((nt)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;



%% Get primary
fid = fopen('flat_ghost.bin','r');
data = fread(fid,[nt n2],'single');
fclose(fid);
fid = fopen('flat_primary_160.bin','r');
origdata = fread(fid,[nt n2],'single');
fclose(fid);
withghost = data(600:800,40);
orig = origdata(600:800,40);


nt = 201;
n2 = 100;
dt = 0.00066666;
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


zmax = 30;
zmin = 30;
tmaxcheck = 1.5*zmax*2/vwater;
fmaxcheck = fmax/4;
fmincheck = fix(1/tmaxcheck);
tmincheck1 = 1./fmaxcheck;
tmincheck2 = 2*zmin/vwater;

tmincheck = 2*zmin/vwater;
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
 tmincheck
 tmaxcheck
 tau2 = (ind -1)*dt
 if ( (tau2>=tmincheck)&&(tau2<=tmaxcheck) )
   tau2
 end
   minustshift = exp(1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2) + epsilon;
   fwithghost = fft(withghost);
   fP =( fwithghost.*F2)./denominator;
   for i = nt:-1:floor(nt/2)
    fP(i) = conj(fP(nt - i + 2));
   end
   
   P = real(ifft(fP));
 
   G = withghost - P;
   

 
 %% Calculate P 
%  weight = abs(P+G);
%  maxweight = max(weight);
%  normweight = weight/maxweight;
%  P = P.*normweight;
 
 
 fP = fft(P);
 absfP = abs(fP);
%  
%  
%  %% scale energy
%  max1 = max(absfwithghost);
%  max2 = max(absfP);
%  scale = max1/max2;
%  fP = fP*scale;
%  absfP = abs(fP);
% 
%  
%  P = real(ifft(fP));
 
 
 P = (addback*withghost + P)/(1.0+addback);
 
 
 
 meancheck = mean(abs(check));
 maxcheck = max(abs(check));
 noisecriteria = maxcheck/meancheck
 
 
 figure(6)
 plot(t,withghost,'blue',t,P,'red',t,orig,'green');
 legend('withghost','P','orig');
 
 
 figure(7);
 plot(tcheck,testmean,'red',tcheck,maxvalue,'green',tcheck,check,'blue');
 legend('testmean','maxvalue','check');
 
 
 figure(8);
 plot(t,absfwithghost,'red',t,absfP,'blue');  
 legend('absfwithghost','absfP');











