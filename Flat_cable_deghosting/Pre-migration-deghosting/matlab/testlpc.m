close all
clear all
clc

iter_cg = 10000;
epsilon = 0.5;

nt = 2001;
n2 = 1001;
dt = 0.004;
tmax = dt*(nt-1);
df = 1./tmax;
fmax = df*(nt-1);

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;

fmin=0;
fmax=100;

fminindex = floor(fmin/df) + 1;
fmaxindex = floor(fmax/df) + 1;

z = 50.0;
tau = 2.0*z/1500;
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;

fid = fopen('primary.bin','r');
synthetic = fread(fid,[nt n2],'float');
fclose(fid);


%% Get primary
primary_sig = synthetic(1:nt,floor(n2/2));
noise = randn(2001,1)/10.0;
primary_sig = primary_sig + 1./cosh(100*(t-1)) +  1./cosh(100*(t-2)) + 1./cosh(100*(t-4)) + 1./cosh(100*(t-5)) + 1./cosh(100*(t-6));
figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F);% + noise;
figure(2)
plot(t,withghost);


%% Frequency domain
fwithghost = fft(withghost);
figure(3)
plot(f,abs(fwithghost));


