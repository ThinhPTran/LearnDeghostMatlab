clear all
close all
clc


%% Deghost parameters
d_src = 16;
d_rec = 32;
vwater = 1500;
refcoe = -0.9; 
recover_eps = 0.02;



%% Input data parameters
nx = 401;
nt = 500;
ntau = nt;
np = 1001;

dt = 0.004;
dx = 8.0;
dtau = dt;
dp = 2.2/(vwater*(np-1));
fp = -1.1/vwater;
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



%% input
fid = fopen('data_with_ghost','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


exinput = zeros(exnt, exnx);
exinput(1:nt,npad+1:npad+nx) = input;


finput = fft(input);
fexinput = fft(exinput);
input_fk = fft2(input);
exinput_fk = fft2(exinput);


fid = fopen('fdata_with_ghost','wb');
fwrite(fid,abs(finput),'single');
fclose(fid);

fid = fopen('data_with_ghost_fk','wb');
fwrite(fid,abs(input_fk),'single');
fclose(fid);


numerator = zeros(size(exinput));
denominator = zeros(size(exinput));
filter1 = zeros(size(exinput));
filter2 = zeros(size(exinput));



%% Generate Filter 
for i = 1:floor(exnf/2)+1
  for j = 1:exnkx
     tmp_value = (exomega(i)*exomega(i))/(vwater*vwater) - exkx(j)*exkx(j);
     if (tmp_value>0)
        numerator(i,j) = 1 + refcoe*exp(1i*2*d_rec*sqrt(tmp_value));
        denominator(i,j) = 1 + refcoe*refcoe + 2*refcoe*cos(2*d_rec*sqrt(tmp_value)) + recover_eps;
        filter1(i,j) = numerator(i,j)/denominator(i,j);
     else
        filter1(i,j) = exp(1i*sqrt(tmp_value));
     end
  end
end

for i = 1:floor(exnf/2)+1
  for j = 1:exnkx
     tmp_value = (exomega(i)*exomega(i))/(vwater*vwater) - exkx(j)*exkx(j);
     if (tmp_value>0)
        numerator(i,j) = 1 + refcoe*exp(1i*2*d_rec*sqrt(tmp_value));
        denominator(i,j) = 1 + refcoe*refcoe + 2*refcoe*cos(2*d_rec*sqrt(tmp_value)) + recover_eps;
        filter2(i,j) = numerator(i,j)/denominator(i,j);
     else
        filter2(i,j) = exp(1i*sqrt(tmp_value));
     end
  end
end



%% Filter 
exdeghost_fk = exinput_fk.*filter1;
exdeghost_fk = exdeghost_fk.*filter2;


% Keep symmetry
for i = 2:floor(exnf/2)+1
  exdeghost_fk(exnf+2-i,1) = conj(exdeghost_fk(i,1));
end


for i = 2:floor(exnf/2)+1
  for j = 2:exnx
    exdeghost_fk(exnt+2-i,exnx+2-j) = conj(exdeghost_fk(i,j));
  end
end

exdeghost = real(ifft2(exdeghost_fk));
exdeghost = cleanup(1.1,1500,exnt,ext,exnf,exf,exomega,exnx,exx,exnp,exxp,exnkx,exkx,exdeghost);
exdeghost_fk = fft2(exdeghost);


deghost = exdeghost(1:nt,npad+1:npad+nx);



%% Output
fid = fopen('exdeghost_fk.bin','wb');
fwrite(fid,abs(exdeghost_fk),'single');
fclose(fid);

fid = fopen('deghost.bin','wb');
fwrite(fid,deghost,'single');
fclose(fid);

fid = fopen('exdeghost.bin','wb');
fwrite(fid,exdeghost,'single');
fclose(fid);



finput=fft(input);
fdeghost=fft(deghost); 

sumfinput=sum(finput,2);
sumfdeghost=sum(fdeghost,2); 

abssumfinput=mag2db(abs(sumfinput)); 
abssumfdeghost=mag2db(abs(sumfdeghost)); 


figure(); 
plot(f,abssumfinput,f,abssumfdeghost); 
legend('abssumfinput','abssumfdeghost'); 








