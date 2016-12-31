clear all
close all
clc



%% Control parameters
vwater = 1500;
recover_eps = 0.01;
p_control = 1.0;
ctr_n = 1.0;
iter_cg = 1000;
res = 5.0;



%% Input data parameters
nx = 101;
nt = 500;
ntau = nt;
np = 1001;

dt = 0.004;
dx = 8.0;
dtau = dt;
dp = 2.0/(vwater*(np+1));
fp = -1.0/vwater+dp;
fx = 100;


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



%% Input data parameters
nx3 = 3;
nt3 = 500;
ntau3 = nt3;
np3 = 201;

dt3 = 0.004;
dx3 = 8.0;
dtau3 = dt3;
dp3 = 2.0/(vwater*(np3+1));
fp3 = -1.0/vwater+dp3;
fx3 = 0;


nf3 = nt3;
nkx3 = nx3;
df3 = 1.0/((nt3)*dt3);
dkx3 = 2.0*pi/((nx3)*dx3);


t3=(0:dt3:(nt3-1)*dt3)';
f3=(0:df3:(nf3-1)*df3)';
kx3=dkx3*[0:floor(nkx3/2) -floor(nkx3/2):-1]';
x3=(fx3:dx3:fx3+(nx3-1)*dx3)';
tau3=(0:dtau3:(ntau3-1)*dtau3)';
p3 = (fp3:dp3:-fp3)';
omega3 = 2.0*pi*f3;



%% Extend on the boundaries
npad3 = 0;
exnx3 = nx3 + 2*npad3;
exnt3 = 2*nt3;
exntau3 = exnt3;
exnp3 = np3;

exdt3 = dt3;
exdx3 = dx3;
exdtau3 = dtau3;
exdp3 = dp3;
exfp3 = fp3;
exfx3 = fx3-npad3*exdx3;

exnf3 = exnt3;
exnkx3 = exnx3;
exdf3 = 1/((exnt3)*exdt3);
exdkx3 = 2.0*pi/((exnx3)*dx3);

ext3=(0:exdt3:(exnt3-1)*exdt3)';
exf3=(0:exdf3:(exnf3-1)*exdf3)';
exkx3=exdkx3*[0:floor(exnkx3/2) -floor(exnkx3/2):-1]';
exx3=(exfx3:exdx3:exfx3+(exnx3-1)*exdx3)';
extau3=(0:exdtau3:(exntau3-1)*exdtau3)';
exxp3 = (exfp3:exdp3:-exfp3)';
exomega3 = 2.0*pi*exf3;




%% Deghost parameters
x0=-20;
x_noff=100;
alpha=atan(1.0/8.0);

d=(x+x0)*tan(alpha);

fmin=1;
fmax=120;
nfb=4;



%% input data
fid = fopen('slant_with_ghost_case1.bin','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


suminput=sum(input,2);
finput=abs(fft(suminput));


figure();
plot(f,finput);
legend('finput');


figure();
imagesc(input);
title('input');


output=zeros(size(input));



%% Forward & Deghost & Backward

tic
for i_iter=2:nx-1
  i_iter
  % Take three traces
  threetrcs=input(:,i_iter-1:i_iter+1);
  threetrcs_taup=taup_fwd(nf3,f3,nx3,x3,np3,p3,threetrcs);
  threetrcs_taup_deghost=taup_deghost(p_control,vwater,d(i_iter), ...
      recover_eps,nt3,dt3,t3,omega3,np3,dp3,p3,threetrcs_taup);
%   threetrcs_taup_deghost_clean=cleanup(threetrcs_taup,...
%       threetrcs_taup_deghost,np3,nt3,dt3,nfb,fmin,fmax,res);
  threetrcs_inv=taup_bwd(iter_cg,nt3,t3,nf3,f3,nx3,x3,np3,p3,threetrcs_taup_deghost);
  output(:,i_iter)=threetrcs_inv(:,2);
end
toc


output1=clean_non_physical_fk(ctr_n,vwater,nt,t,nf,f,omega,nx,x,np,p,nkx,kx,output);


output_fk=fft2(output);
output1_fk=fft2(output1);


% figure();
% plot(t,threetrcs(:,2));
% legend('threetrcs_input');
% 
% 
% figure();
% imagesc(threetrcs_taup);
% title('threetrcs\_taup');
% 
% 
% figure();
% imagesc(threetrcs_taup_deghost);
% title('threetrcs\_taup\_deghost');


% figure();
% imagesc(threetrcs_taup_deghost_clean);
% title('threetrcs\_taup\_deghost\_clean');



%% Check output
sumoutput=sum(output,2);
foutput=abs(fft(sumoutput));


figure();
imagesc(output);
title('output');


figure();
imagesc(output1);
title('output1');


figure();
imagesc(abs(output_fk));
title('output_fk');


figure();
imagesc(abs(output1_fk));
title('output1_fk');


figure();
plot(f,foutput);
legend('foutput');


%% Output data to bin
fid=fopen('deghost.bin','wb');
fwrite(fid,output,'single');
fclose(fid);


fid=fopen('deghost_fk.bin','wb');
fwrite(fid,abs(output1_fk),'single');
fclose(fid);


in=threetrcs(:,2);
out=threetrcs_inv(:,2);

figure();
plot(t,in,t,out);
legend('in','out');


