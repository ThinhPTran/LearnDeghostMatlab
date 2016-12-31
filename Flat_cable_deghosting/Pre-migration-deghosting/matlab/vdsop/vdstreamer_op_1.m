close all
clear all
clc


designednt = 1001;
designeddt = 0.004;
designeddf = 1./((designednt-1)*designeddt);
designedtmax = (designednt-1)*designeddt;
designedfmax = (designednt-1)*designeddf;

designedt = (0:designeddt:designedtmax)';
designedf = (0:designeddf:designedfmax)';
designedomega = 2*pi*designedf;


designedz = 9;
designedtau = 2.0*designedz/1500;
designedamplitude = 100;
designedindex = floor(designedtau/designeddt) + floor(designednt/2) + 1;

tmpF = zeros(designednt, 1);
tmpF(floor(designednt/2) + 1) = 1.0;
tmpF(designedindex) = -1.0;


%% Get primary
noise = designedamplitude*randn(designednt,1)*0.05;
wlet = ricker(40,0.004);
primary_sig =  zeros(designednt,1) +designedamplitude*1./cosh(100*(designedt-0.3));
for i = 1:13
   primary_sig(i+10) = 1.0*designedamplitude*wlet(i); 
   primary_sig(i+150) = 0.8*designedamplitude*wlet(i);
   primary_sig(i+300) = 1.3*designedamplitude*wlet(i);
   primary_sig(i+450) = 0.5*designedamplitude*wlet(i);
   primary_sig(i+600) = 1.1*designedamplitude*wlet(i);
   primary_sig(i+750) = 0.3*designedamplitude*wlet(i);
   primary_sig(i+900) = 1.0*designedamplitude*wlet(i);
end

figure(1)
plot(designedt,primary_sig);


addback = 0.1;

ntwindow = designednt;
dtwindow = designeddt;
dfwindow = 1./((ntwindow-1)*dtwindow);
tmaxwindow = (ntwindow-1)*dtwindow;
fmaxwindow = (ntwindow-1)*dfwindow;

twindow = (0:dtwindow:tmaxwindow)';
fwindow = (0:dfwindow:fmaxwindow)';
omegawindow = 2*pi*fwindow;


zmax = 75;
zmin = 15;
tmaxcheck = 1.5*zmax*2/1500;
fmaxcheck = fmaxwindow/4;
fmincheck = fix(1/tmaxcheck);
tmincheck1 = 1./fmaxcheck;
tmincheck2 = 2*zmin/1500;

tmincheck = max([tmincheck1 tmincheck2]);
tmaxcheck;



denominatore = zeros(ntwindow,1);
F = zeros(ntwindow,1);
withghost = zeros(ntwindow,1);
fwithghost = zeros(ntwindow,1);
absfwithghost = zeros(ntwindow,1);
fP = zeros(ntwindow,1);
P = zeros(ntwindow,1);
G = zeros(ntwindow,1);


%% Synthetize data with ghost
tmpsig = convolution2(primary_sig,tmpF)  + noise;
withghost = tmpsig;
figure(2)
plot(designedt,withghost);


tic


fwithghost = fft(withghost);
absfwithghost = abs(fwithghost);


figure(3)
plot(designedf,absfwithghost);


inittau = 0.0;
currenttau=inittau;
maxiter = 100000;
epsilon = 0.0000001;
dt = designeddt;


r = model(withghost,omegawindow,dt,currenttau);
err_new = r'*r


for i_iter = 1:5
    
    i_iter

    J = (model(withghost,omegawindow,dt,currenttau+dt) - model(withghost,omegawindow,dt,currenttau))/dt;

    mean(r)
    mean(J)
    
    
    delta = -mean(r)/mean(J)

    currenttau = currenttau + delta
    
    
    r = model(withghost,omegawindow,dt,currenttau);
    err_old = err_new;
    err_new = r'*r;
    
    rms = sqrt(err_new);
    frac = abs(err_new - err_old)/abs(err_old);
    
    if rms < epsilon ; break; end;
    if frac < epsilon ; break; end;
    
end

currenttau



store = zeros(ntwindow,1);
% 
% for i_iter = 1:10
%     tmptau = (i_iter-1)*dt
%     tmpstore = model(withghost,omegawindow,dt,tmptau);
%     store(i_iter) = mean(tmpstore);
% end
% 
% 
%     figure(4)
%     plot(designedt,store);






toc


