clear all
close all
clc

eps = 0.0001;
iter_cg = 1000;
vwater = 1500;


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
x=(fx:dx:-fx)';
tau=(0:dtau:(ntau-1)*dtau)';
p = (fp:dp:-fp)';
omega = 2.0*pi*f;


%% Extend on the boundaries
npad = 50;
exnx = nx + 2*npad;
exnt = nt;
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
exx=(exfx:exdx:-exfx)';
extau=(0:exdtau:(exntau-1)*exdtau)';
exxp = (exfp:exdp:-exfp)';
exomega = 2.0*pi*exf;



fid = fopen('receiver_data_with_ghost_n1_500_taup.bin','r');
input = fread(fid,[exnt exnp],'single');
fclose(fid);



finput = fft(input);


fid = fopen('receiver_data_with_ghost_n1_500_taup_deghost.bin','r');
input1 = fread(fid,[exnt exnp],'single');
fclose(fid);



finput1 = fft(input1);





A = zeros(np,nx);
output = zeros(nt,nx);
foutput = zeros(nf,nx);
exoutput = zeros(nt,exnx);
exfoutput = zeros(nf,exnx);
output1 = zeros(nt,nx);
foutput1 = zeros(nf,nx);
exoutput1 = zeros(nt,exnx);
exfoutput1 = zeros(nf,exnx);
tmp_in = zeros(np,1);
tmp_out = zeros(exnx,1);


tic

%% Inverse for original data
for i_iter = 1:floor(nf/2)
  
  tmp_in = transpose(finput(i_iter,:));
    
  for j_iter = 1:exnp
    for k_iter = 1:exnx
      A(j_iter,k_iter) = exp(1i*2.0*pi*exf(i_iter)*exx(k_iter)*exxp(j_iter)); 
    end
  end
  
  tmp_out = conj_grad_solve(A,tmp_in,iter_cg);

   
  exfoutput(i_iter,:) = transpose(tmp_out);
  
end


foutput = exfoutput(:,npad+1:npad+nx);


for i_iter = 2:floor(nf/2)+1
  foutput(nf-i_iter+2,:) = conj(foutput(i_iter,:)); 
end


output = real(ifft(foutput));


toc


fid = fopen('receiver_data_with_ghost_n1_500_invtaup.bin','w');
fwrite(fid,output,'single');
fclose(fid);


fid = fopen('data_with_ghost','rb');
input = fread(fid,[nt nx],'single');
fclose(fid);


err = sqrt(mean2((input - output).^2))


tic

%% Inverse for deghosted data
for i_iter = 1:floor(nf/2)
  
  tmp_in = transpose(finput1(i_iter,:));
    
  for j_iter = 1:exnp
    for k_iter = 1:exnx
      A(j_iter,k_iter) = exp(1i*2.0*pi*exf(i_iter)*exx(k_iter)*exxp(j_iter)); 
    end
  end
  
  tmp_out = conj_grad_solve(A,tmp_in,iter_cg);

   
  exfoutput1(i_iter,:) = transpose(tmp_out);
  
end

foutput1 = exfoutput1(:,npad+1:npad+nx);

for i_iter = 2:floor(nf/2)+1
  foutput1(nf-i_iter+2,:) = conj(foutput1(i_iter,:));
  exfoutput1(nf-i_iter+2,:) = conj(exfoutput1(i_iter,:));
end


output1 = real(ifft(foutput1));
exoutput1 = real(ifft(exfoutput1));

foutput1 = fft(output1);
output1_fk = fft2(output1);
exoutput1_fk = fft2(exoutput1);

toc


%% Filter on FK domain

% Transform into FK domain

fid = fopen('receiver_data_with_ghost_n1_500_invtaup_fk_before_filtering.bin','wb');
fwrite(fid,abs(output1_fk),'single');
fclose(fid);


fid = fopen('exreceiver_data_with_ghost_n1_500_invtaup_fk_before_filtering.bin','wb');
fwrite(fid,abs(exoutput1_fk),'single');
fclose(fid);


% Filtering
disp('Filtering...');
for i_iter = 1:floor(exnf/2)+1
  for j_iter = 1:exnkx
     if ((exomega(i_iter)*exomega(i_iter))/(vwater*vwater) - exkx(j_iter)*exkx(j_iter)<=0)
       exoutput1_fk(i_iter,j_iter) = 0;
     end
  end
end


% Keep the symmetry
disp('Keep the symmetry...');
for i_iter = 2:floor(exnf/2)+1
  exoutput1_fk(exnf+2-i_iter,1) = conj(exoutput1_fk(i_iter,1));
end


for i_iter = 2:floor(exnf/2)+1
  for j_iter = 2:exnkx
    exoutput1_fk(exnf+2-i_iter,exnkx+2-j_iter) = conj(exoutput1_fk(i_iter,j_iter));
  end
end


% Output to check FK domain
fid = fopen('exreceiver_data_with_ghost_n1_500_invtaup_fk_after_filtering.bin','wb');
fwrite(fid,abs(exoutput1_fk),'single');
fclose(fid);



% Recalculate results
exoutput1 = real(ifft2(exoutput1_fk));
output1 = exoutput1(:,npad+1:npad+nx);
foutput1 = fft(output1);
output1_fk = fft2(output1);


% for i_iter = 1:nt
%   for j_iter = 1:nx
%     if ((t(i_iter) < 6.0e-4*x(j_iter) + 0.18) || (t(i_iter) < -6.0e-4*x(j_iter) + 0.18))
%       output1(i_iter, j_iter) = 0.0;
%     end
%   end
% end
% 
% output1(1:108,:) = 0.0;



%% Output the final results
fid = fopen('exreceiver_data_with_ghost_n1_500_invtaup_deghost.bin','w');
fwrite(fid,exoutput1,'single');
fclose(fid);


fid = fopen('receiver_data_with_ghost_n1_500_invtaup_deghost.bin','w');
fwrite(fid,output1,'single');
fclose(fid);


fid = fopen('freceiver_data_with_ghost_n1_500_invtaup_deghost.bin','w');
fwrite(fid,abs(foutput1),'single');
fclose(fid);


fid = fopen('receiver_data_with_ghost_n1_500_invtaup_fk_deghost.bin','w');
fwrite(fid,abs(output1_fk),'single');
fclose(fid);



