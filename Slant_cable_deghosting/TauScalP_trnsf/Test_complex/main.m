clear all
close all
clc



%% Control parameters
vwater=1500;
recover_eps=0.15;
p_control=1.0;
ctr_n=1.1;
iter_cg=1000;
res=5.0;



%% Input data parameters
nx=240;
nt=501;
ntau=nt;
np=1501;

dt=0.004;
dx=12.5;
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



%% Acquision parameters
% x0=-20;
% x_noff=100;
alpha=atan(0.9375/12.5);

% d=(x+x0)*tan(alpha);



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



%% input data
fid=fopen('slantcable_complex_09375.bin','r');
input=fread(fid,[nt nx],'single');
fclose(fid);


figure();
imagesc(input);
title('input');


input_tausclp=tausclp_fwd(nf,f,nx,x,np,p,input);


figure();
imagesc(input_tausclp);
title('input\_tausclp');


tic
input_tausclp_dg=taup_slant_dg(recover_eps,alpha,vwater,nt,dt,npprime,dpprime,input_tausclp);
toc


figure();
imagesc(input_tausclp);
title('input\_tausclp');


figure();
imagesc(input_tausclp_dg);
title('input\_tausclp\_dg');


tic
input_inv=tausclp_bwd(iter_cg,nt,t,nf,f,nx,x,np,p,input_tausclp_dg);
input_inv1=tausclp_bwd(iter_cg,nt,t,nf,f,nx,x,np,p,input_tausclp);
toc


output=clean_non_physical_fk(ctr_n,vwater,nt,t,nf,f,omega,nx,x,np,p,nkx,kx,input_inv);


output_fk=abs(fft2(output));


diff=input-input_inv1;


%% View and write to file
fid=fopen('deghost.bin','wb');
fwrite(fid,output,'single');
fclose(fid);

fid=fopen('input_inv.bin','wb'); 
fwrite(fid,input_inv1,'single'); 
fclose(fid); 

fid=fopen('diff.bin','wb'); 
fwrite(fid,diff,'single'); 
fclose(fid);  



figure();
imagesc(output);
title('output');


figure();
imagesc(output_fk);
title('output\_fk');


figure(); 
imagesc(input_inv1); 
title('input inv'); 





















