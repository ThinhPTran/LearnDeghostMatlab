close all
clear all
clc

iter_cg = 10000;
epsilon = 0.5;

nt = 2001;
n2 = 1001;
dt = 0.004;
tmax = (nt-1)*dt;
df = 1./tmax;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';


z = 50.0;
tau = 2.0*z/1500;
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;

% fid = fopen('primary.bin','r');
% synthetic = fread(fid,[nt n2],'float');
% fclose(fid);


%%Get primary
primary_sig = synthetic(1:nt,floor(n2/2));
noise = randn(nt,1)/10.0;
primary_sig = primary_sig + 1./cosh(100*(t-2)) + 1./cosh(100*(t-4));
figure(1)
plot(t,primary_sig);
%%%%%%%%%%%%%


%%Synthetize data with ghost
withghost = convolution2(primary_sig,F);% + noise;
figure(2)
plot(t,withghost);
%%%%%%%%%%%%%


%% Calculate spectrum of the synthesized data
fwithghost = fft(withghost);
figure(3)
plot(f,abs(fwithghost));
%%%%%%%%%%%%%


%% Find the primary
b = zeros(nt,1);

fprimary = find_fprimary(fwithghost,b,iter_cg,epsilon);













