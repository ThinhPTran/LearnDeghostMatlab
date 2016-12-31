%% Test analyse_spec_fwd

clear all
close all
clc



%% Control parameters
vwater = 1500;
recover_eps = 0.05;
p_control = 1.0;
ctr_n = 1.0;
iter_cg = 1000;
res = 5.0;



%% Input data parameters
nx = 401;
nt = 500;
ntau = nt;
np = 1001;

dt = 0.004;
dx = 8.0;
dtau = dt;
dp = 2.0/(vwater*(np-1));
fp = -1.0/vwater;
fx = -1600;


nf = nt;
nkx = nx;
df = 1.0/((nt)*dt);
dkx = 2.0*pi/((nx)*dx);


t=(0:dt:(nt-1)*dt)';
f=(0:df:(nf-1)*df)';
kx=dkx*[0:floor(nkx/2) -floor(nkx/2):-1]';
x=(fx:dx:fx+(nx-1)*dx)';
tau=(0:dtau:(ntau-1)*dtau)';
p = (fp:dp:-fp)';
omega = 2.0*pi*f;



%% Extend on the boundaries
npad = 100;
exnx = nx + 2*npad;
exnt = 2*nt;
exntau = exnt;
exnp = np;

exdt = dt;
exdx = dx;
exdtau = dtau;
exdp = dp;
exfp = fp;
exfx = fx-npad*exdx;

exnf = exnt;
exnkx = exnx;
exdf = 1/((exnt)*exdt);
exdkx = 2.0*pi/((exnx)*dx);

ext=(0:exdt:(exnt-1)*exdt)';
exf=(0:exdf:(exnf-1)*exdf)';
exkx=exdkx*[0:floor(exnkx/2) -floor(exnkx/2):-1]';
exx=(exfx:exdx:exfx+(exnx-1)*exdx)';
extau=(0:exdtau:(exntau-1)*exdtau)';
exxp = (exfp:exdp:-exfp)';
exomega = 2.0*pi*exf;



%% Deghost parameters
d_src = 16.0;
d_rec = 32.0;

fmin=1
fmax=125
nfb=15



%% input data
fid=fopen('data_with_ghost','r');
input=fread(fid,[nt nx],'single');
fclose(fid);

fid=fopen('output1.bin','r');
output=fread(fid,[nt nx],'single');
fclose(fid);



%% Take one trace
in_trc=input(:,201);
out_trc=output(:,201);


figure();
plot(in_trc);
legend('in\_trc');


figure();
plot(out_trc);
legend('out\_trc');


[tspecin,tmp]=analyse_spec_fwd1(in_trc,1,nt,dt,nfb,fmin,fmax,res);
[tspecout,tmp]=analyse_spec_fwd1(out_trc,1,nt,dt,nfb,fmin,fmax,res);


for i_iter = 1:nfb
%   figure();
%    plot(tspecin(:,i_iter));
%   plot(f,abs(tspecin(:,i_iter)),f,abs(tspecout(:,i_iter)));
%   legend('band in','band out');
end


for i_iter = 1:nfb
  maxin=max(abs(tspecin(:,i_iter)));
  maxout=max(abs(tspecout(:,i_iter)));
  ratio=maxout/maxin;
  if (ratio>0)
    tspecout(:,i_iter)=tspecout(:,i_iter)/ratio; 
  end
end


tmp_tspecin=zeros(size(tspecin));
tmp_tspecout=zeros(size(tspecout));


for i_iter = 1:nfb
  tmp_tspecin(:,i_iter)=real(ifft(tspecin(:,i_iter)));
  tmp_tspecout(:,i_iter)=real(ifft(tspecout(:,i_iter)));
end
  

in_trc_inv=sum(tmp_tspecin,2);
out_trc_inv=sum(tmp_tspecout,2);


figure();
plot(t,in_trc,t,in_trc_inv);
legend('in\_trc','in\_trc\_inv');


figure();
plot(t,out_trc,t,out_trc_inv);
legend('out\_trc','out\_trc\_inv');



%% Compare input and output spec
fin_trc=abs(fft(in_trc));
fin_trc_inv=abs(fft(in_trc_inv));
fout_trc=abs(fft(out_trc));
fout_trc_inv=abs(fft(out_trc_inv));



figure();
plot(f,fin_trc,f,fin_trc_inv);
legend('fin\_trc','fin\_trc\_inv');


figure();
plot(f,fin_trc,f,fout_trc,f,fout_trc_inv);
legend('fin\_trc','fout\_trc','fout\_trc\_inv');


figure();
plot(t,in_trc,t,out_trc,t,out_trc_inv);
legend('in\_trc','out\_trc','out\_trc\_inv');









