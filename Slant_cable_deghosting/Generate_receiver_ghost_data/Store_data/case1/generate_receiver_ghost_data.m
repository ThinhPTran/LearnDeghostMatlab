close all
clear all
clc


%% Control parameters
vwater=1500;
recover_eps=0.15;
p_control=1.0;
ctr_n=1.1;
iter_cg=1000;
res=5.0;



%% Input data parameters
nx=101;
nt=500;
ntau=nt;
np=1501;

dt=0.004;
dx=8.0;
dtau=dt;
dp=2.0/(vwater*(np+1));
fp=-1.0/vwater+dp;
fx=100;


nf=nt;
nkx=nx;
df=1.0/((nt)*dt);
dkx=2.0*pi/((nx)*dx);


t=(0:dt:(nt-1)*dt)';
f=(0:df:(nf-1)*df)';
kx=dkx*[0:floor(nkx/2) -floor(nkx/2):-1]';
x=(fx:dx:fx+(nx-1)*dx)';
tau=(0:dtau:(ntau-1)*dtau)';
p=(fp:dp:-fp)';
omega=2.0*pi*f;



fid=fopen('slant_realcable.bin','rb');
realcable=fread(fid,[nt nx],'single');
fclose(fid);


fid=fopen('slant_mirrcable.bin','rb');
mirrcable=fread(fid,[nt nx],'single');
fclose(fid);


fid=fopen('flat_realcable_case1.bin','rb');
flatcable=fread(fid,[nt, nx],'single');
fclose(fid);



%% FB mute
for i_iter = 1:nt
  for j_iter = 1:nx
    if ((i_iter < j_iter + 90))
      realcable(i_iter, j_iter) = 0.0;
      mirrcable(i_iter,j_iter) = 0.0;
    end
  end
end



data_with_ghost=realcable-mirrcable;



%% Generate flat
alpha=atan(1.05/8.0);



%% Deghost parameters

% Calculate pprime
pmin=p(1)
pmax=p(np)
pprimemin=pmin*vwater/sqrt(1-pmin*pmin*vwater*vwater)
pprimemax=pmax*vwater/sqrt(1-pmax*pmax*vwater*vwater)

pprimemin=-4;
pprimemax=4;

% ptest=1/vwater
% pprimetest=ptest/sqrt(1-ptest*ptest*vwater*vwater)

npprime=np;
dpprime=(pprimemax-pprimemin)/(npprime-1);
pprime=(pprimemin:dpprime:pprimemax)';

for i_iter=1:np
  p(i_iter)=pprime(i_iter)/(vwater*sqrt(1+pprime(i_iter)*pprime(i_iter))); 
end



figure();
plot(1:np,vwater*p);
title('p');


realcable_tausclp=tausclp_fwd(nf,f,nx,x,np,p,realcable);


figure();
imagesc(realcable_tausclp);
title('real\_cable\_tausclp');


tic
realcablerot=taup_slant_rot(recover_eps,alpha,vwater,nt,dt,npprime,dpprime,realcable_tausclp);
toc


figure();
imagesc(realcablerot);
title('real\_cable\_tausclp\_rot');


tic
flatcable_tmp=tausclp_bwd(iter_cg,nt,t,nf,f,nx,x,np,p,realcablerot);
toc


output=clean_non_physical_fk(ctr_n,vwater,nt,t,nf,f,omega,nx,x,np,p,nkx,kx,flatcable_tmp);


output_fk=abs(fft2(output));


err=output-flatcable;



%% View data
figure();
imagesc(realcable);
title('realcable');

figure()
imagesc(mirrcable);
title('mirrcable');

figure();
imagesc(data_with_ghost);
title('data\_with\_ghost');

figure();
imagesc(flatcable_tmp);
title('flat\_cable\_tmp');

figure();
imagesc(output);
title('output');

figure();
imagesc(flatcable);
title('flat\_cable');

figure();
imagesc(err);
title('err');



%% Write output
fid=fopen('slant_with_ghost_case1.bin','wb');
fwrite(fid,data_with_ghost,'single');
fclose(fid);


fid=fopen('flat_cable_noghost_case1.bin','wb');
fwrite(fid,output,'single');
fclose(fid);







