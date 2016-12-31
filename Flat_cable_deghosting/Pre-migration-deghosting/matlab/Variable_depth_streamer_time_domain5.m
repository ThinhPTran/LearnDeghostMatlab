close all
clear all
clc

iter_cg = 10000
epsilon = 0.0001

nt = 2001;
n2 = 1001;
dt = 0.004
tmax = dt*(nt-1)
df = 1./tmax
fmax = df*(nt-1)

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;

fmin=0;
fmax=100;

fminindex = floor(fmin/df) + 1;
fmaxindex = floor(fmax/df) + 1;

z = 50.0;
tau = 2.0*z/1500
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;

fid = fopen('primary.bin','r');
synthetic = fread(fid,[nt n2],'float');
fclose(fid);


%% Get primary
primary_sig = synthetic(1:nt,floor(n2/2));
noise = randn(2001,1)/100.0;
primary_sig = primary_sig + 1./cosh(100*(t-1)) +  1./cosh(100*(t-2)) + 1./cosh(100*(t-4)) + 1./cosh(100*(t-5)) + 1./cosh(100*(t-6));
figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F) + noise;
figure(2)
plot(t,withghost);


%% Initial value for optimization process
 tau2 = tau


 %% Calculate P second way.
 minustshift = exp(1i*tau2*omega);
 a = 1.0;
 F2 = 1 - a*minustshift;
 denominator = 1 + a*a -  2*a*cos(omega*tau2) + 0.1;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P2 = real(ifft(fP));
 figure(3)
 plot(t,P2);

 
 %% Calculate G
 G2 = withghost - P2;
 figure(4)
 plot(t,G2);
 
 
 %% Find tau
 fP2 = fft(P2);
 fG2 = fft(G2);
 
 
 ftau = find_tau(fG2,-fP2,iter_cg,epsilon);
 
 for i = nt:-1:floor(nt/2)
  ftau(i) = conj(ftau(nt - i + 2));
 end
 
 tau = real(ifft(ftau));
 
 figure(5)
 plot(f,abs(ftau));


