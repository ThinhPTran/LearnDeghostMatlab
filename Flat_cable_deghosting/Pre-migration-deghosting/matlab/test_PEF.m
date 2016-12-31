close all
clear all
clc

iter_cg = 10000;
epsilon = 0.01;

nt = 2001;
n2 = 1001;
dt = 0.004;

t = (0:dt:8)';

z = 50.0;
tau = 2.0*z/1500;
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -0.8;

fid = fopen('primary.bin','r');
synthetic = fread(fid,[nt n2],'float');
fclose(fid);


%%Get primary
primary_sig = synthetic(1:nt,floor(n2/2));
noise = randn(nt,1)/10.0;
primary_sig = primary_sig + 1./cosh(100*(t-2)) + 1./cosh(100*(t-4));
figure(1)
plot(t,primary_sig);
%%%%%%%%%%%%%


%%Synthetize data with ghost
withghost = convolution2(primary_sig,F) + noise;
figure(2)
plot(t,withghost);
%%%%%%%%%%%%%


tmp = forward_PEF(withghost,iter_cg,epsilon);
forwardfilter = zeros(nt,1);
forwardfilter(floor(nt/2) + 1) = tmp(1);
forwardfilter(floor(nt/2) + 2) = tmp(2);

forwardpredictwithghost = convolution2(withghost,forwardfilter);
forwardpredictwithghost = circshift(forwardpredictwithghost,[1 0]);

tmp = backward_PEF(withghost,iter_cg,epsilon);
backwardfilter = zeros(nt,1);
backwardfilter(floor(nt/2) + 1) = tmp(1);
backwardfilter(floor(nt/2) + 2) = tmp(2);

backwardpredictwithghost = convolution2(withghost,backwardfilter);
backwardpredictwithghost = circshift(backwardpredictwithghost,[-2 0]);

predictwithghost = (forwardpredictwithghost + backwardpredictwithghost)/2.0;

figure(3);
plot(t,predictwithghost);





