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
% nx = 101;
% nt = 500;
nx=240;
nt=501;
ntau = nt;
np = 1001;

dt = 0.004;
dx = 8.0;
dtau = dt;
dp = 2.0/(vwater*(np-1));
fp = -1.0/vwater;
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



%% Extend on the boundaries
npad = 0;
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
x0=-20;
x_noff=100;
alpha=atan(1.0/8.0);

d_src = 16.0;

fmin=1;
fmax=120;
nfb=4;



%% input data
fid = fopen('slantcable_09375.bin','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


input_fk = fft2(input);

fid = fopen('data_with_ghost_fk','wb');
fwrite(fid,abs(input_fk),'single');
fclose(fid);


exinput = zeros(exnt, exnx);
exinput(1:nt,npad+1:npad+nx)=input;


%% View data
figure();
imagesc(input);
title('input');

figure()
imagesc(abs(input_fk));
title('taup\_fk');



%% TauP transform
disp('Slant TauScP transform');
tic
exinput_slanttaup=slant_tauscp_fwd(vwater,x0,x_noff,alpha,exnf,exf,exnx,exx,exnp,exxp,exinput);
toc

fid = fopen('exreceiver_data_with_ghost_n1_1000_taup.bin','wb');
fwrite(fid,exinput_slanttaup,'single');
fclose(fid);

%% View data
figure();
imagesc(exinput_slanttaup);
title('exinput\_slanttaup');



%% Get extmp_inv
disp('Inverse on TauP domain for extmp_inv0');
tic
extmp_inv0=slant_tauscp_bwd(iter_cg,vwater,x0,x_noff,alpha,exnf,exf,exnx,exx,exnp,exxp,exinput_slanttaup);
toc

fid = fopen('exreceiver_data_with_ghost_n1_1000_invtaup.bin','wb');
fwrite(fid,extmp_inv0,'single');
fclose(fid);

%% View data
figure();
imagesc(extmp_inv0);
title('extmp\_inv0');


%% Get missing part after TauP transform
disp('Get missing part after TauP transform');
missing_part=exinput-extmp_inv0;

fid = fopen('exreceiver_data_with_ghost_n1_1000_missingpart.bin','wb');
fwrite(fid,missing_part,'single');
fclose(fid);

%% View data
figure();
imagesc(missing_part);
title('missing\_part');


%% Deghost in TauP domain
disp('Deghost on TauP domain');
tic
exinput_slanttaup_deghost=slant_taup_deghost(p_control,vwater,x0,x_noff,alpha,recover_eps,exnt,exdt,ext,exomega,exnp,exdp,exxp,exinput_slanttaup);
toc

exerr_slanttaup = abs(exinput_slanttaup_deghost - exinput_slanttaup);

%% View data
figure();
imagesc(exinput_slanttaup_deghost);
title('exinput\_slanttaup\_deghost');

figure();
imagesc(exerr_slanttaup);
title('exerr\_slanttaup');


%% Clean up
disp('Cleanup on TauP domain');
tic
exinput_slanttaup_deghost_clean=cleanup(exinput_slanttaup,exinput_slanttaup_deghost,exnp,exnt,exdt,nfb,fmin,fmax,res);
toc


exerr_slanttaup_clean = abs(exinput_slanttaup_deghost_clean - exinput_slanttaup);

fid = fopen('exerr_slanttaup_clean.bin','wb');
fwrite(fid,exerr_slanttaup_clean,'single');
fclose(fid);


fid = fopen('receiver_data_with_ghost_n1_1000_taup_deghost.bin','wb');
fwrite(fid,exinput_slanttaup_deghost,'single');
fclose(fid);

fid = fopen('receiver_data_with_ghost_n1_1000_taup_deghost_clean.bin','wb');
fwrite(fid,exinput_slanttaup_deghost_clean,'single');
fclose(fid);


exinput_slanttaup_deghost_sum=sum(exinput_slanttaup_deghost,2);
exinput_slanttaup_deghost_clean_sum=sum(exinput_slanttaup_deghost_clean,2);
fexinput_slanttaup_deghost_sum=abs(fft(exinput_slanttaup_deghost_sum));
fexinput_slanttaup_deghost_clean_sum=abs(fft(exinput_slanttaup_deghost_clean_sum));


figure();
plot(exf,fexinput_slanttaup_deghost_sum,exf,fexinput_slanttaup_deghost_clean_sum);
legend('fexinput\_slanttaup\_deghost\_sum','fexinput\_slanttaup\_deghost\_clean\_sum');



%% Inverse TauP transform
disp('Inverse on TauP domain');
tic 
extmp_inv=slant_tauscp_bwd(iter_cg,vwater,x0,x_noff,alpha,exnf,exf,exnx,exx,exnp,exxp,exinput_slanttaup_deghost_clean);
toc

tmp_exoutput=clean_non_physical_fk(ctr_n,vwater,exnt,ext,exnf,exf,exomega,exnx,exx,exnp,exxp,exnkx,exkx,extmp_inv);


%% View data
figure();
imagesc(tmp_exoutput);
title('tmp_exoutput');


%% Scale back for fix scale display
disp('Scale back for fix scale compare');
tic
max1=max(max(abs(input)));
max2=max(max(abs(tmp_exoutput)));
ratio=max2/max1;
tmp_exoutput1=tmp_exoutput/ratio;
toc

output1=tmp_exoutput1(1:nt,npad+1:npad+nx);
output1_fk=fft2(output1);


%% Add back missing part
disp('Add back missing part');
tic
tmp_exoutput2=tmp_exoutput1+missing_part;
toc

output2=tmp_exoutput2(1:nt,npad+1:npad+nx);
output2_fk=fft2(output2);


% disp('Reshape the spectrum');
% tic
% output3=reshapespec(input,output2,nx,nt,dt,nfb,fmin,fmax,res);
% toc
% output3_fk = fft2(output3);


fid = fopen('output1.bin','wb');
fwrite(fid,output1,'single');
fclose(fid);


fid = fopen('output1_fk.bin','wb');
fwrite(fid,abs(output1_fk),'single');
fclose(fid);


fid = fopen('exoutput1.bin','wb');
fwrite(fid,tmp_exoutput,'single');
fclose(fid);


fid = fopen('output2.bin','wb');
fwrite(fid,output2,'single');
fclose(fid);


fid = fopen('output2_fk.bin','wb');
fwrite(fid,abs(output2_fk),'single');
fclose(fid);


% fid = fopen('output3.bin','wb');
% fwrite(fid,output3,'single');
% fclose(fid);
% 
% 
% fid = fopen('output3_fk.bin','wb');
% fwrite(fid,abs(output3_fk),'single');
% fclose(fid);


sum_input=sum(input,2);
sum_output1=sum(output1,2);
sum_output2=sum(output2,2);
% sum_output3=sum(output3,2);


fsum_input=abs(fft(sum_input));
fsum_output1=abs(fft(sum_output1));
fsum_output2=abs(fft(sum_output2));
% fsum_output3=abs(fft(sum_output3));


% figure();
% plot(f,fsum_input,f,fsum_output1,f,fsum_output2,f,fsum_output3);
% legend('fsum\_input','fsum\_output1','fsum\_output2','fsum\_output3');
% 
% 
% 
% mean2(abs(output3))









