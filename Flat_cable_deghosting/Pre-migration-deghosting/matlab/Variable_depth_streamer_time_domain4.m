close all
clear all
clc

iter_cg = 10000;
epsilon = 0.001;

nt = 101
n2 = 101
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
primary_sig =  1./cosh(100*(t-0.1)) + 1./cosh(100*(t-0.3));% + noise;
figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F) + noise;
figure(2)
plot(t,withghost);


%% Create Mirror data
mirror_data = -withghost;
figure(3)
plot(t,mirror_data);



%% Fourier transform
fwithghost = fft(withghost);
figure(4)
plot(f,abs(fwithghost));

absfwithghost=abs(fwithghost);


 %% Calculate forward filter
 tau2 = 0.0267
 eps = 0.2
 tshift = exp(-1i*tau2*omega);
 a = 1.0;
 F1 = 1 - a*tshift;

 for i = nt:-1:floor(nt/2)
  F1(i) = conj(F1(nt - i + 2));
 end

 for i = 1:nt
  if ( abs(F1(i)) < eps )
    F1(i) = eps;
  end
 end
 
 
figure(5)
plot(f,abs(F1));


 %% Find P
 for i = 1:nt
   fP(i) = fwithghost(i)/F1(i); 
 end


 P = real(ifft(fP));
 
 figure(6)
 plot(t,P);
 
 figure(7)
 plot(f,abs(fP));

 
