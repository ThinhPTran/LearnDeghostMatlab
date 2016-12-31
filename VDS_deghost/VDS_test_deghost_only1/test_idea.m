close all
clear all
clc

%% Basic params
eps = 0.1;
vwater = 1500;
amplitude = 1;


nt = 2001;
n2 = 648;
dt = 0.004;

df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


%  Clean-up params
fmin_hz=0.01;
fmax_hz=(nt-2)*df;
res=4.0;
nfb=35;
ntwindow=30;
wshift=10;
thres=0.5;



%% t for checking
zrec = zeros(n2,1);
zsrc = zeros(n2,1); 
zout = zeros(n2,1); 
tau = zeros(n2,1);


% read rec_elev
zrec=load('vds_data1_rec_elev.txt','-ascii'); 
zsrc=load('vds_data1_sou_elev.txt','-ascii'); 


tau=2.0*zrec(1)/vwater



%% Get data
noise = 0.0001*amplitude*randn(nt,1);


fid = fopen('vds_data1_3.bin','rb');
datain = fread(fid,[nt n2],'single');
fclose(fid);



P = zeros(size(datain));


testiter=50


tic

%% Deghost 
for iter = testiter:testiter
iter
P(:,iter) = deghostfunc(datain(:,iter), zrec(iter), zsrc(iter), vwater, eps, nt, dt); 

%zrecout;
%zscrout; 


%% Do CWT
disp('Do CWT'); 


[time_spectral, freq_hz] = cwt_morlet_fwd(P(:,iter),1,nt,dt,nfb,fmin_hz,fmax_hz,res);



%% Balance
disp('Do balance'); 
time_spectral_out = tspecbalance(time_spectral,nfb,nt,dt,ntwindow,wshift,thres);



%% Calculate output
output=zeros(nt,1); 

for j_iter=1:nfb
    output(:,1)=output(:,1)+real(time_spectral_out(:,j_iter));
end


% % keep max
maxinput = mean(abs(datain(:,iter))); 
maxoutput = mean(abs(output)); 

if (maxoutput>0) 
  output  = output*maxinput/maxoutput; 
end


P(:,iter) = output; 

end


 
figure();
plot(t,datain(:,testiter),'blue',t,P(:,testiter),'red');
legend('withghost','Primary recovered');



finput=abs(fft(datain(:,testiter)));
foutput=abs(fft(P(:,testiter)));

finputdb = mag2db(finput);
foutputdb = mag2db(foutput);


figure();
plot(f,finput,f,foutput);
legend('finput','foutput'); 
title('Raw spectrum'); 


figure();
plot(f,finputdb,f,foutputdb); 
legend('finputdb','foutputdb');
title('Spectrum DB'); 


fid = fopen('data_in.bin','wb');
fwrite(fid,datain,'single');
fclose(fid);


fid = fopen('data_out.bin','wb');
fwrite(fid,P,'single');
fclose(fid);


% figure();
% imagesc(datain); 
% 
% figure();
% imagesc(P); 











