clear all
close all
clc



%% Input parameters
nx = 401;
nkx = 401;
nt = 500;
nf = 500;


dt = 0.004;
dx = 8;
domega = 2.0*pi/((nt)*dt);
df = 1.0/(nt*dt);
dkx = 2.0*pi/((nx)*dx);


t=(0:dt:(nt-1)*dt)';
f=df*[0:floor(nf/2) -floor(nf/2)+1:-1]';
omega = 2*pi*f;
x=(0:dx:(nx-1)*dx)';
kx = dkx*[0:floor(nkx/2) -floor(nkx/2):-1];



%% Deghost parameters
d_rec = 32;
vwater = 1500;
recover_eps = 0.01;



%% input
fid = fopen('data_with_ghost','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


finput = fft(input);
input_fk = fft2(input);


fid = fopen('fdata_with_ghost','wb');
fwrite(fid,abs(finput),'single');
fclose(fid);

fid = fopen('data_with_ghost_fk','wb');
fwrite(fid,abs(input_fk),'single');
fclose(fid);


numerator = zeros(size(input));
denominator = zeros(size(input));
filter = zeros(size(input));



%% Generate Filter 
for i = 1:floor(nf/2)+1
  for j = 1:nkx
     tmp_value = (omega(i)*omega(i))/(vwater*vwater) - kx(j)*kx(j);
     tmp_value1 = (omega(i)*omega(i))/(6000*6000) - kx(j)*kx(j);
     if (tmp_value>0)
       if (tmp_value1>0) 
         filter(i,j) = 0; 
       else
         numerator(i,j) = 1 - exp(1i*2*d_rec*sqrt(tmp_value));
         denominator(i,j) = 2 - 2*cos(2*d_rec*sqrt(tmp_value)) + recover_eps;
         filter(i,j) = numerator(i,j)/denominator(i,j);
       end
     else
        filter(i,j) = exp(1i*sqrt(tmp_value));
     end
  end
end



%% Filter 
deghost_fk = input_fk.*filter;


% Keep symmetry
for i = 2:floor(nf/2)+1
  deghost_fk(nf+2-i,1) = conj(deghost_fk(i,1));
end


for i = 2:floor(nf/2)+1
  for j = 2:nx
    deghost_fk(nt+2-i,nx+2-j) = conj(deghost_fk(i,j));
  end
end

deghost = real(ifft2(deghost_fk));



%% Output
fid = fopen('deghost_fk.bin','wb');
fwrite(fid,abs(deghost_fk),'single');
fclose(fid);


fid = fopen('deghost.bin','wb');
fwrite(fid,deghost,'single');
fclose(fid);





