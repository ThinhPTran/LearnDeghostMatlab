clear all
close all
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
omega = 2*pi*f;

z = 50.0;
tau = 2.0*z/1500
index = floor(tau/dt) + floor(nt/2) + 1;

tau = floor(tau/dt)*dt

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;

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
withghost = convolution2(primary_sig,F);% + noise;
figure(2)
plot(t,withghost);
%%%%%%%%%%%%%


%%Create Mirror data
shift_back = zeros(nt,1);
index1 = -floor(tau/dt) + floor(nt/2) + 1;
shift_back(index1) = 1.0;
mirror_data = convolution2(shift_back,withghost);
figure(3)
plot(t,mirror_data);
%%%%%%%%%%%%%



%%Find the Primary (Deghosting) 

tic;

   
 F1 = cal_filter(mirror_data, withghost, iter_cg, epsilon);
% F1 = cal_filter(withghost, withghost, iter_cg, epsilon);


figure(4)
plot(t,F1);

%%%%%%%%%%%%%


%% Create ghost filter

maxind = 1;
maxvalue = abs(F1(1));
for i = 2:nt
    if ( abs(F1(i)) > maxvalue ) 
       maxindex = i;
       maxvalue = abs(F1(i));
    end
end

F1 = zeros(nt,1);
F1(floor(nt/2)+1) = 1.0;
F1(maxindex) = -1.0;


%% Calculate P

% P = cal_primary(F1, withghost, iter_cg, epsilon);
%  
%  
%  figure(5);
%  plot(t,P);
 
 
 %% Calculate P second way.
%  tau2 = (maxindex - floor(nt/2) - 1)*dt
%  minustshift = exp(1i*tau2*omega);
%  F2 = 1 - minustshift;
%  denominator = 2 -  2*cos(omega*tau2) + 1e-15;
%  fwithghost = fft(withghost);
%  fP =( fwithghost.*F2)./denominator;
%  for i = nt:-1:floor(nt/2)
%   fP(i) = conj(fP( - i + 2));
%  end
% 
%  P2 = real(ifft(fP));
%  figure(5)
%  plot(t,P2);
 
 
%  figure(6)
%  plot(f,abs(fP));
 
 
  %% Calculate P third way.
 tau2 = (maxindex - floor(nt/2) - 1)*dt
 f = (0:0.125:250)';
 omega = 2*pi*f;
 shift = exp(-1i*tau2*omega);
 F2 = 1 - shift;
 fwithghost = fft(withghost);
 fP =( fwithghost./(F2 + 0.1) );
 for i = 2001:-1:floor(2001/2)
  fP(i) = conj(fP(2001 - i + 2));
 end

 P2 = real(ifft(fP));
 figure(5)
 plot(t,P2);
 

toc;



