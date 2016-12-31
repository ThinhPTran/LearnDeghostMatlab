close all
clear all
clc

iter_cg = 10000;
epsilon = 0.1;
stable_p = 1.0;
addback = 0.1;
weight = 0;

nt= 4000;
n2 = 576;
dt = 0.002;
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;



%% Get primary
fid = fopen('filter9mghost.bin','r');
data = fread(fid,[nt n2],'single');
fclose(fid);

dataout = zeros(nt,n2);


for tridx = 1:n2
    
   withghost = data(:,tridx);
    
   tau2 = 0.012;

   minustshift = exp(1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2);
   stable_p = 1.0*mean(abs(denominator));
   denominator = denominator + stable_p;
   fwithghost = fft(withghost);
   fP =( fwithghost.*F2)./denominator;
   for i = nt:-1:floor(nt/2)
    fP(i) = conj(fP(nt - i + 2));
   end

   P2 = real(ifft(fP));
 
   G = withghost - P2;
   

   dataout(:,tridx) = P2;
    
end

fid = fopen('output.bin','w');
fwrite(fid,dataout,'single');
fclose(fid);

