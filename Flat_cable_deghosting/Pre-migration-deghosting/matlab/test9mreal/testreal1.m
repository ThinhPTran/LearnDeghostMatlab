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
check = zeros(ntcheck,1);

findcheck = fix(tmincheck/dt) + 1;
lindcheck = fix(tmaxcheck/dt) + 1;


weight = wshift./ntwindow;
denominatorstore = zeros(ntwindow,lindcheck);
F2store = zeros(ntwindow,lindcheck);
denominatore = zeros(ntwindow,1);
F2 = zeros(ntwindow,1);
withghost = zeros(ntwindow,1);
fwithghost = zeros(ntwindow,1);
absfwithghost = zeros(ntwindow,1);
fP = zeros(ntwindow,1);
P2 = zeros(ntwindow,1);
G2 = zeros(ntwindow,1);



for i = 1:lindcheck
    tmptau = (i-1)*dt;
    denominatorstore(:,i) = 2 -  2*cos(omega*tmptau) + 2.0;
    minustshift = exp(1i*tmptau*omega);
    F2store(:,i) = 1 - minustshift;
end





for tridx = 1:n2
    for twidx = 1:nw
        
tridx
twidx
        
fidx = (twidx-1)*wshift + 1;
lidx = (twidx-1)*wshift + ntwindow;

withghost = data(fidx:lidx,tridx);


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);

average = mean(abs(withghost));
averagespec = mean(absfwithghost);


for iter = 1:lindcheck
 F2 = F2store(:,iter);
 denominator = denominatorstore(:,iter);
 fP =( fwithghost.*F2)./denominator;
 for i = ntwindow:-1:floor(ntwindow/2)
  fP(i) = conj(fP(ntwindow - i + 2));
 end

 P2 = real(ifft(fP));
 
 maxvalue(iter) = max(P2)/average;

end


for iter = 1:lindcheck
 F2 = F2store(:,iter);
 denominator = denominatorstore(:,iter);
 fP =( fwithghost.*F2)./denominator;
 for i = ntwindow:-1:floor(ntwindow/2)
  fP(i) = conj(fP(ntwindow - i + 2));
 end

 P2 = real(ifft(fP));
 G2 = withghost - P2;
 
 test = G2 + circshift(P2,[iter-1 0]);
 testmean(iter) = 1./mean(abs(test)/average);

end


for i = 1:ntcheck
    check(i) = 1.5*testmean(i)+maxvalue(i);
end


[tmp ind] = max(check);

%% Recover
 tau2 = (ind -1)*dt;
 if ( abs(tau2 - 0.012) < 0.04) 
   F2 = F2store(:,ind);
   denominator = 2 -  2*cos(omega*tau2) + 0.2;
   fP =( fwithghost.*F2)./denominator;
   for i = ntwindow:-1:floor(ntwindow/2)
    fP(i) = conj(fP(ntwindow - i + 2));
   end

   P2 = real(ifft(fP));
   
 else 
   P2 = withghost;
 end
 
 
 taustore(twidx,tridx) = tau2;
 
 P2 = (addback*withghost + P2)/(1.0+addback);
 

 
 for i = 1:ntwindow
     dataout(fidx + i - 1,tridx) = dataout(fidx + i - 1,tridx) + weight*P2(i);
 end
 
 
    end
end


fid = fopen('output.bin','w');
fwrite(fid,dataout,'single');
fclose(fid);