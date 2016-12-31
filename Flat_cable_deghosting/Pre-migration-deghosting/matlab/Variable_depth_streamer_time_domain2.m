%close all
%clear all
%clc

iter_cg = 10000;
epsilon = 0.5;

nt = 2001;
n2 = 1001;
dt = 0.004;
df = 0.125;

t = (0:dt:8)';
f = (0:df:250)';
omega = 2*pi*f;

fmin=0;
fmax=100;

fminindex = floor(fmin/df) + 1;
fmaxindex = floor(fmax/df) + 1;


z = 20.0;
tau = 2.0*z/1500
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;

fid = fopen('primary.bin','r');
synthetic = fread(fid,[nt n2],'float');
fclose(fid);


%%Get primary
primary_sig = synthetic(1:nt,floor(n2/2));
noise = randn(2001,1)/10.0;
primary_sig = primary_sig + 1./cosh(100*(t-2)) + 1./cosh(100*(t-4));% + noise;
figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F);% + noise;
figure(2)
plot(t,withghost);


%%% Denoise
%fwithghost = fft(withghost);
%
%for i=fmaxindex+1:nt
  %fwithghost(i)=0;
%end
%
%for i = nt:-1:floor(nt/2)
  %fwithghost(i) = conj(fwithghost(nt - i + 2));
%end
%
%withghost = real(ifft(fwithghost));
%
%figure(3);
%plot(t,withghost);


%% Create input for crosscorr2
corrwithghost(1) = withghost(i);
corrwithghost(nt) = withghost(nt);
for i = 3:nt
  a = (withghost(i-1) - withghost(i-2))/dt;
  err = abs(withghost(i-2) + a*2.0*dt  - withghost(i));
  corrwithghost(i) = withghost(i)*withghost(i)./err;
end

figure(4)
plot(t,corrwithghost);


%% Create Mirror data
corrmirror_data = -corrwithghost;
figure(5)
plot(t,corrmirror_data);



%% Find the Primary (Deghosting) 

tic;

   
 F1 = crosscorr2(mirror_data,withghost);


figure(6)
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

figure(7);
plot(t,F1);


%% Calculate P

% P = cal_primary(F1, withghost, iter_cg, epsilon);
%  
%  
%  figure(6);
%  plot(t,P);
 
 
 %% Calculate P second way.
 tau2 = (maxindex - floor(nt/2) - 1)*dt
 a = 1.0;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - a*minustshift;
 denominator = 2 -  2*cos(omega*tau2) + 1e-10;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P2 = real(ifft(fP));
 figure(8)
 plot(t,P2);
 

toc;