close all
clear all
clc

iter_cg = 10000;
epsilon = 0.1;
vwater = 1500;
amplitude = 1;


nt = 101;
n2 = 101;
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
primary_sig =  zeros(nt,1);
for i = 1:13
   primary_sig(i+10) = 1.0*amplitude*wlet(i); 
end

figure(1)
plot(t,primary_sig);
legend('primary');



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F)  + noise;
%withghost = noise;
figure(2)
plot(t,withghost);
legend('withghost');



%% Calculate fft
absfprimary_sig = abs(fft(primary_sig));
anglefprimary_sig = angle(fft(primary_sig));
absfwithghost = abs(fft(withghost));
anglefwithghost = angle(fft(withghost));


frecovered_sig = absfprimary_sig.*exp(1i*anglefwithghost);
recovered_sig = real(ifft(frecovered_sig));


figure(3)
plot(f,absfprimary_sig,f,absfwithghost);
legend('fprimary','fwithghost');


figure(4)
plot(t,primary_sig,t,withghost,t,recovered_sig);
legend('primary\_sig','withghost','recovered\_sig');









