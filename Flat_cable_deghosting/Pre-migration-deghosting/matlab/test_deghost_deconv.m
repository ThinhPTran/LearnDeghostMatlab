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
index = floor(tau/dt) + 1;

ntcheck_min = floor(taumin_check/dt)+1;
ntcheck_max = floor(taumax_check/dt)+1;

F = zeros(nt, 1);
F(1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = 0.1*amplitude*randn(nt,1);
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



withghost1=conv(primary_sig,F);
withghost2=withghost1+0.05*randn(size(withghost1));
withghost2(1000:end)=0.0;


figure(3)
plot(withghost2);


deghost_deconv=deconv(withghost2,F);


figure(4)
plot(deghost_deconv);




