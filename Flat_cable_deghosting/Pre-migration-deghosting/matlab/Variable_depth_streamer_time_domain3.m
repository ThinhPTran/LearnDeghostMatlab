close all
clear all
clc

iter_cg = 10000;
epsilon = 0.1;

nt = 51
n2 = 1001
dt = 0.004
df = 1./((nt-1)*dt)
tmax = (nt-1)*dt
fmax = (nt-1)*df

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


z = 20.0;
tau = 2.0*z/1500
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = randn(nt,1)/20.0;
primary_sig =  1./cosh(100*(t-0.07));% + 1./cosh(100*(t-4));% + noise;
figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F);% + noise;
figure(2)
plot(t,withghost);


%% Create Mirror data
mirror_data = -withghost;
figure(3)
plot(t,mirror_data);




%% Find the Primary (Deghosting) 

tic;

   
 F1 = crosscorr2(withghost, -withghost);


figure(4)
plot(t,F1);



%% Create ghost filter

maxind = 1;
maxvalue = F1(1);
for i = 2:nt
    if ( F1(i) > maxvalue ) 
       maxindex = i;
       maxvalue = F1(i);
    end
end

if (maxindex < (floor(nt/2) + 1)) 
  tmp = floor(nt/2) + 1 - maxindex;
  maxindex = floor(nt/2) + 1 + tmp;
end

F1 = zeros(nt,1);
F1(floor(nt/2)+1) = 1.0;
F1(maxindex) = -1.0;

figure(5);
plot(t,F1);
 

 
%% Calculate P

P = cal_primary(F1, withghost, iter_cg, epsilon);
 
 
figure(6);
plot(t,P);


%% Calculate G

G = withghost - P;

figure(7);
plot(t,G);
 

toc;

 
