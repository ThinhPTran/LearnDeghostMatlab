close all
clear all
clc

taustore = zeros(1000,1);

for k=1:1


nt = 51;
dt = 0.004;
df = 1./((nt)*dt);
tmax = (nt)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;



z = 6;
tau = 2.0*z/1500;
%tau = 0.02
amplitude = 100;
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = amplitude*randn(nt,1)*0.05;
wlet = ricker(40,0.004);
primary_sig1 =  zeros(nt,1);
primary_sig2 =  zeros(nt,1);
for i = 1:13
  primary_sig1(i+15) = 1.0*amplitude*wlet(i); 
%   primary_sig2(i+19) = 1.0*amplitude*wlet(i);
end

primary_sig = primary_sig1 + primary_sig2;
%  primary_sig = 1./cosh(100*(t-0.06)) - 1./cosh(100*(t-0.06-tau));

% figure(1)
% plot(t,primary_sig);


addback = 0.0;
epsilon = 0.00000001;
ntwindow = nt;
dtwindow = dt;
dfwindow = 1./((ntwindow)*dtwindow);
tmaxwindow = (ntwindow-1)*dtwindow;
fmaxwindow = (ntwindow-1)*dfwindow;

twindow = (0:dtwindow:tmaxwindow)';
fwindow = (0:dfwindow:fmaxwindow)';
omegawindow = 2*pi*fwindow;


zmax = 50;
zmin = 6;
tmaxcheck = zmax*2.0/1500;
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
tmp = convolution2(primary_sig,F)+noise;
withghost = tmp;
% figure(2)
% plot(t,withghost);


tic


for i = 1:ntcheck
    tmptau = (i-1)*dtwindow;
    denominatorstore(:,i) = 2 -  2*cos(omegawindow*tmptau) + 2.0;
    minustshift = exp(1i*tmptau*omegawindow);
    F2store(:,i) = 1 - minustshift;
end


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);
% figure(3)
% plot(f,absfwithghost);


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
    check(i) = testmean(i)+maxvalue(i);
end


[tmp ind] = max(check);

%% Recover
 tau2 = (ind -1)*dtwindow; 
 
 
 
 %% Calculate P, 1st way
 F2 = F2store(:,ind);
 denominator = 2 -  2*cos(omegawindow*tau2) + epsilon;
 fP =( fwithghost.*F2)./denominator;
 for i = ntwindow:-1:floor(ntwindow/2)
   fP(i) = conj(fP(ntwindow - i + 2));
 end

 
 
 P = real(ifft(fP));
 G = withghost-P;
 
 %% Calculate P 
%  weight = abs(P+G);
%  maxweight = max(weight);
%  normweight = weight/maxweight;
%  P = P.*normweight;
 
 
 fP = fft(P);
 absfP = abs(fP);
 
 
%  max1 = max(absfwithghost);
%  max2 = max(absfP);
%  scale = max1/max2;
%  
%  fP = fP*scale;
%  absfP = abs(fP);
% 
%  
%  P = real(ifft(fP));
 
 
 P = (addback*withghost + P)/(1.0+addback);
 
 
 toc
 
 
 tau
 tau2
 
 
 taustore(k)=tau2;
 
 
 meancheck = mean(abs(check));
 maxcheck = max(abs(check));
 noisecriteria = maxcheck/meancheck;
 
 
 figure(4)
 plot(twindow,withghost,'blue',twindow,P,'red',twindow,primary_sig,'green',twindow,G,'yellow');
 hleg = legend('input data','deghosted data','primary','ghost');
 set(hleg,'FontAngle','italic','TextColor',[.3 .2 .1])
 saveas(gcf,'30mnoise.eps','eps2c');
 
 
 figure(5);
 plot(tcheck,testmean,'red',tcheck,maxvalue,'green',tcheck,check,'blue');
 
 
 figure(6);
 plot(twindow,absfwithghost,'red',twindow,absfP,'blue');
 
 
end
 
% figure(7)
% x = (0.02:0.0001:0.04);
% [n xout] =  hist(taustore,x);
% normn = n./sum(n);
% bar(xout,normn,'histc');
% xlabel('$tau$','interpreter','latex','fontsize',18);
% ylabel('$Percentage$','interpreter','latex','fontsize',18);
% title('Normalized histogram for found tau in two objectives case Namp/Sigamp = 0.05');
% set(gca,'fontsize',14); % increase font size
% saveas(gcf,'30mnoise_statistic.eps','eps2c');