close all
clear all
clc

iter_cg = 10000;
epsilon = 0.1;
stable_p = 0.1

nt = 51
n2 = 1001
dt = 0.004
df = 1./((nt-1)*dt)
tmax = (nt-1)*dt
fmax = (nt-1)*df

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


z = 25;
tau = 2.0*z/1500
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = randn(nt,1)/15.0;
wlet = ricker(40,0.004);
primary_sig =  zeros(nt,1) +1./cosh(100*(t-0.14));
for i = 1:13
   primary_sig(i+10) = wlet(i); 
end

figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F)  + noise;
figure(2)
plot(t,withghost);


maxvalue = zeros(nt,1);

for iter = 1:floor(nt/2)+1

 tau2 = (iter-1)*dt;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2) + stable_p;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 P2 = real(ifft(fP));
 
 maxvalue(iter) = max(P2);

end


[tmp ind] = max(maxvalue);
tau2 = (ind-1)*dt;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2) + stable_p;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = nt:-1:floor(nt/2)
  fP(i) = conj(fP(nt - i + 2));
 end

 
 P2 = real(ifft(fP));
 
 maxP2 = max(P2);
 minP2 = min(P2);
 
 if ( abs(maxP2) < 1.5*abs(minP2)) 
    maxvalue(ind) = 0;
    
    [max ind] = max(maxvalue);
    tau2 = (ind-1)*dt;
    minustshift = exp(1i*tau2*omega);
    F2 = 1 - minustshift;
    denominator = 2 -  2*cos(omega*tau2) + stable_p;
    fwithghost = fft(withghost);
    fP =( fwithghost.*F2)./denominator;
    for i = nt:-1:floor(nt/2)
    fP(i) = conj(fP(nt - i + 2));
    end
    P2 = real(ifft(fP));
 end
 
 
 G = withghost - P2;
 
 
 tau2


figure(3)
plot(t,withghost,'blue',t,P2,'red',t,G,'green');






 
 









