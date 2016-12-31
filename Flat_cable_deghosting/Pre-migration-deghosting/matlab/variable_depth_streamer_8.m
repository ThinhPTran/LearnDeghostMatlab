close all
clear all
clc

iter_cg = 10000;
epsilon = 0.1;
stable_p = 1.0

nt = 101
n2 = 1001
dt = 0.004
df = 1./((nt-1)*dt)
tmax = (nt-1)*dt
fmax = (nt-1)*df

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


z = 50;
tau = 2.0*z/1500
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = randn(nt,1)/15.0;
wlet = ricker(40,0.004);
primary_sig =  zeros(nt,1) +1./cosh(100*(t-0.3));
for i = 1:13
   primary_sig(i+10) = 1.0*wlet(i); 
%    primary_sig(i+150) = 1.0*wlet(i);
%    primary_sig(i+300) = 1.0*wlet(i);
%    primary_sig(i+450) = 1.0*wlet(i);
%    primary_sig(i+600) = 1.0*wlet(i);
%    primary_sig(i+750) = 1.0*wlet(i);
%    primary_sig(i+900) = 1.0*wlet(i);
end

figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F)  + noise;
figure(2)
plot(t,withghost);


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);
figure(3)
plot(f,absfwithghost);


zmax = 50
zmin = 15;
tmaxcheck = 1.5*zmax*2/1500 
fmaxcheck = fmax/4
tmincheck1 = 1./fmaxcheck;
tmincheck2 = 2*zmin/1500;
tmincheck = max([tmincheck1 tmincheck2])

ntcheck = fix(tmaxcheck/dt) + 1;
tcheck = 0:dt:(ntcheck-1)*dt;
maxvalue = zeros(ntcheck,1);
baseonf = zeros(ntcheck,1);
check = zeros(ntcheck,1);

findcheck = fix(tmincheck/dt) + 1;
lindcheck = fix(tmaxcheck/dt) + 1;


for iter = findcheck:lindcheck

 tau2 = (iter-1)*dt;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2) + stable_p;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P2 = real(ifft(fP));
 
 maxvalue(iter) = max(P2);

end


figure(4);
plot(tcheck,maxvalue);


for i = 2:ntcheck
   tmpt = (i-1)*dt;
   tmpf = 1./t;
   find = fix((1./tmpt)/df) + 1;
   baseonf(i) = (1./(absfwithghost(find) + 0.000001));
end


figure(5);
plot(tcheck,baseonf);


for i = 1:ntcheck
    check(i) = maxvalue(i)*baseonf(i);
end



figure(6);
plot(tcheck,check,'blue',tcheck,maxvalue,'red',tcheck,baseonf,'green');


[tmp ind] = max(check);

%% Recover
 tau2 = (ind -1)*dt
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2) + stable_p;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P2 = real(ifft(fP));
 
 G = withghost - P2;
 
figure(7)
plot(t,withghost,'blue',t,P2,'red',t,G,'green');
 
 
 



