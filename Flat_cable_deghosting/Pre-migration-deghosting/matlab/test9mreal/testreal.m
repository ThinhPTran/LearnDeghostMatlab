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


ntwindow = 50;
wshift = 5;
nw = fix((nt - ntwindow)/wshift) + 1;
n2 = 576;
dt = 0.002;
df = 1./((ntwindow-1)*dt);
tmax = (ntwindow-1)*dt;
fmax = (ntwindow-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


dataout = zeros(nt,n2);
taustore = zeros(nw,n2);


for tridx = 1:n2
    for twidx = 1:nw
        
tridx
twidx
        
fidx = (twidx-1)*wshift + 1;
lidx = (twidx-1)*wshift + ntwindow;

withghost = data(fidx:lidx,tridx);


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);


zmax = 50;
zmin = 15;
tmaxcheck = 1.5*zmax*2/1500;
fmaxcheck = fmax/4;
fmincheck = fix(1/tmaxcheck);
tmincheck1 = 1./fmaxcheck;
tmincheck2 = 2*zmin/1500;

tmincheck = max([tmincheck1 tmincheck2]);
tmaxcheck;
fmincheck;
fmaxcheck;

ntcheck = fix(tmaxcheck/dt) + 1;
tcheck = 0:dt:(ntcheck-1)*dt;
testmean = zeros(ntcheck,1);
maxvalue = zeros(ntcheck,1);
baseonf = zeros(ntcheck,1);
check = zeros(ntcheck,1);

findcheck = fix(tmincheck/dt) + 1;
lindcheck = fix(tmaxcheck/dt) + 1;


average = mean(abs(withghost));
averagespec = mean(absfwithghost);


for iter = 1:lindcheck

 tau2 = (iter-1)*dt;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2);
 stable_p = 1.0*mean(abs(denominator));
 denominator = denominator + 1.0*stable_p;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = ntwindow:-1:floor(ntwindow/2)
  fP(i) = conj(fP(ntwindow - i + 2));
 end

 P2 = real(ifft(fP));
 
 maxvalue(iter) = max(P2)/average;

end


for iter = 1:lindcheck

 tau2 = (iter-1)*dt;
 minustshift = exp(1i*tau2*omega);
 F2 = 1 - minustshift;
 denominator = 2 -  2*cos(omega*tau2);
 stable_p = 1.0*mean(abs(denominator));
 denominator = denominator + 1.0*stable_p;
 fwithghost = fft(withghost);
 fP =( fwithghost.*F2)./denominator;
 for i = ntwindow:-1:floor(ntwindow/2)
  fP(i) = conj(fP(ntwindow - i + 2));
 end

 P2 = real(ifft(fP));
 G2 = withghost - P2;
 
 test = G2 + circshift(P2,[iter-1 0]);
 testmean(iter) = 1./mean(abs(test)/average);

end


for i = 2:ntcheck
   tmpt = (i-1)*dt;
   tmpf = 1./t;
   find = fix((1./tmpt)/df) + 1;
   baseonf(i) = (1./(absfwithghost(find)/averagespec + 0.000001));
end


for i = 1:ntcheck
    check(i) = 1.5*testmean(i)+maxvalue(i);
end


[tmp ind] = max(check);

%% Recover
 tau2 = (ind -1)*dt;
 if ( abs(tau2 - 0.012) < 0.04) 
   minustshift = exp(1i*tau2*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*tau2);
   stable_p = 0.1*mean(abs(denominator));
   denominator = denominator + stable_p;
   fwithghost = fft(withghost);
   fP =( fwithghost.*F2)./denominator;
   for i = ntwindow:-1:floor(ntwindow/2)
    fP(i) = conj(fP(ntwindow - i + 2));
   end

   P2 = real(ifft(fP));
 
   G = withghost - P2;
   
 else 
   P2 = withghost;
   G = withghost - P2;
 end
 
 
 taustore(twidx,tridx) = tau2;
 
 P2 = (addback*withghost + P2)/(1.0+addback);
 

 
 for i = 1:ntwindow
     weight = wshift./ntwindow;
     dataout(fidx + i - 1,tridx) = dataout(fidx + i - 1,tridx) + weight*P2(i);
 end
 
 
    end
end


fid = fopen('output.bin','w');
fwrite(fid,dataout,'single');
fclose(fid);










