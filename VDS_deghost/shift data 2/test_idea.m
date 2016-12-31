close all
clear all
clc



eps = 2.0;
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



%% t for checking
z = 10
fpeak = 80;


tau = 2.0*z/vwater


%% Get primary
noise = 0.05*amplitude*randn(nt,1);
wlet = rickerwave(fpeak,0.2,t);
primary_sig = zeros(nt,1) +amplitude*1./cosh(800*(t-0.2));


figure;
plot(t,primary_sig);
legend('primary_sig');


%% Synthetize data with ghost
fprimary_sig = fft(primary_sig);
fwithghost = fprimary_sig.*(1-exp(-1i*tau*omega));
for j = nt:-1:floor(nt/2)-1
    fwithghost(j) = conj(fwithghost(nt - j + 2));
end


withghost = real(ifft(fwithghost))+noise;%+0.4*amplitude*1./cosh(800*(t-0.3));


figure
plot(t,withghost);
legend('withghost');


%% Deghost 
[P zout] = deghostfunc(withghost, z, vwater, eps, nt, dt); 


zout

 
figure
plot(t,withghost,'blue',t,P,'red',t,primary_sig,'green');
legend('withghost','Primary recovered','Primary');


















