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
nx = 240;
nt = 501;
ntau = nt;
np = 1001;

dt = 0.004;
dx = 12.5;
dtau = dt;
dp = 2.0/(vwater*(np+1));
fp = -1.0/vwater+dp;
fx = 0;


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
fmax=120
nfb=50



%% input data
fid = fopen('data_with_ghost','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


input_fk = fft2(input);

fid = fopen('data_with_ghost_fk','wb');
fwrite(fid,abs(input_fk),'single');
fclose(fid);

figure();
imagesc(abs(input_fk));
title('input\_fk');


exinput = zeros(exnt, exnx);
exinput(1:nt,npad+1:npad+nx)=input;
 


%% TauP transform
disp('TauP transform');
tic
exinput_taup = taup_fwd(exnf,exf,exnx,exx,exnp,exxp,exinput);
toc

% fid = fopen('exreceiver_data_with_ghost_n1_1000_taup.bin','wb');
% fwrite(fid,exinput_taup,'single');
% fclose(fid);

figure();
imagesc(exinput_taup);
title('exinput\_taup');


%% Get extmp_inv
maxinput=max(max(abs(exinput_taup)));
exinput_taup_tmp = exinput_taup + 0.001*maxinput*randn(size(exinput_taup));
disp('Inverse on TauP domain for extmp_inv0');
tic
extmp_inv0=taup_bwd(iter_cg,exnt,ext,exnf,exf,exnx,exx,exnp,exxp,exinput_taup);
toc

% fid = fopen('exreceiver_data_with_ghost_n1_1000_invtaup.bin','wb');
% fwrite(fid,extmp_inv0,'single');
% fclose(fid);

figure();
imagesc(extmp_inv0);
title('extmp\_inv0');


%% Get missing part after TauP transform
disp('Get missing part after TauP transform');
missing_part=exinput-extmp_inv0;

% fid = fopen('exreceiver_data_with_ghost_n1_1000_missingpart.bin','wb');
% fwrite(fid,missing_part,'single');
% fclose(fid);

figure();
imagesc(missing_part);
title('missing\_part');


%% Deghost in TauP domain
disp('Deghost on TauP domain');
tic
exinput_taup_deghost=taup_deghost(p_control,vwater,d_src,d_rec,recover_eps,exnt,exdt,ext,exomega,exnp,exdp,exxp,exinput_taup);
toc


exerr_taup = abs(exinput_taup_deghost - exinput_taup);

% fid = fopen('exerr_taup.bin','wb');
% fwrite(fid,exerr_taup,'single');
% fclose(fid);

figure();
imagesc(exerr_taup);
title('exerr\_taup after deghost - before deghost in TauP domain');


%% Clean up
disp('Cleanup on TauP domain');
tic
exinput_taup_deghost_clean=cleanup(exinput_taup,exinput_taup_deghost,exnp,exnt,exdt,nfb,fmin,fmax,res);
toc


exerr_taup_clean = abs(exinput_taup_deghost_clean - exinput_taup);


% fid = fopen('exerr_taup_clean.bin','wb');
% fwrite(fid,exerr_taup_clean,'single');
% fclose(fid);

figure();
imagesc(exerr_taup_clean);
title('exerr\_taup\_after clean up - before clean up in TauP domain');


% fid = fopen('receiver_data_with_ghost_n1_500_taup_deghost.bin','wb');
% fwrite(fid,exinput_taup_deghost,'single');
% fclose(fid);

figure();
imagesc(exinput_taup_deghost);
title('exinput\_taup\_deghost');


% fid = fopen('receiver_data_with_ghost_n1_500_taup_deghost_clean.bin','wb');
% fwrite(fid,exinput_taup_deghost_clean,'single');
% fclose(fid);

figure();
imagesc(exinput_taup_deghost_clean);
title('exinput\_taup\_deghost\_clean');


exinput_taup_deghost_sum=sum(exinput_taup_deghost,2);
exinput_taup_deghost_clean_sum=sum(exinput_taup_deghost_clean,2);
fexinput_taup_deghost_sum=abs(fft(exinput_taup_deghost_sum));
fexinput_taup_deghost_clean_sum=abs(fft(exinput_taup_deghost_clean_sum));


figure();
plot(exf,fexinput_taup_deghost_sum,exf,fexinput_taup_deghost_clean_sum);
legend('fexinput\_taup\_deghost\_sum','fexinput\_taup\_deghost\_clean\_sum');



%% Inverse TauP transform
disp('Inverse on TauP domain');
tic
extmp_inv = taup_bwd(iter_cg,exnt,ext,exnf,exf,exnx,exx,exnp,exxp,exinput_taup_deghost_clean);
toc


figure();
imagesc(extmp_inv);
title('extmp\_inv');


tmp_exoutput = clean_non_physical_fk(ctr_n,vwater,exnt,ext,exnf,exf,exomega,exnx,exx,exnp,exxp,exnkx,exkx,extmp_inv);


figure();
imagesc(tmp_exoutput);
title('tmp\_exoutput');


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


disp('Reshape the spectrum');
tic
output3=reshapespec(input,output2,nx,nt,dt,nfb,fmin,fmax,res);
toc
output3_fk=fft2(output3);


disp('Reshape: focus on high frequencies');
tic
output4=reshapespec1(input,output2,nx,nt,dt,nfb,fmin,fmax,res);
toc
output4_fk=fft2(output4);



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


fid = fopen('output3.bin','wb');
fwrite(fid,output3,'single');
fclose(fid);


fid = fopen('output3_fk.bin','wb');
fwrite(fid,abs(output3_fk),'single');
fclose(fid);


sum_input=sum(input,2);
sum_output1=sum(output1,2);
sum_output2=sum(output2,2);
sum_output3=sum(output3,2);
sum_output4=sum(output4,2);


fsum_input=abs(fft(sum_input));
fsum_output1=abs(fft(sum_output1));
fsum_output2=abs(fft(sum_output2));
fsum_output3=abs(fft(sum_output3));
fsum_output4=abs(fft(sum_output4));


fsum_input=mag2db(fsum_input);
fsum_output1=mag2db(fsum_output1);
fsum_output2=mag2db(fsum_output2);
fsum_output3=mag2db(fsum_output3);
fsum_output4=mag2db(fsum_output4);


figure();
plot(f,fsum_input,f,fsum_output1,f,fsum_output2,f,fsum_output3);
legend('fsum\_input','fsum\_output1','fsum\_output2','fsum\_output3');

figure();
plot(f,fsum_input,f,fsum_output1,f,fsum_output3);
legend('fsum\_input','fsum\_output1','fsum\_output3');

figure();
plot(f,fsum_input,f,fsum_output1,f,fsum_output4);
legend('fsum\_input','fsum\_output1','fsum\_output4');



mean2(abs(output3))
mean2(abs(output4))









