clear all
close all
clc


threshold_ratio = 0.3;
vwater = 1500;
recover_eps = 0.001;


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
npad = 0;
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



%% Deghost parameters
pmax = p(np);
d = 32;
zmin_check = 20;
zmax_check = 60;

taumin_check = 2.0*zmin_check/vwater;
taumax_check = 2.0*zmax_check/vwater;

ntcheck_min = floor(taumin_check/dt)+1;
ntcheck_max = floor(taumax_check/dt)+1;




%% input
fid = fopen('data_with_ghost','rb');
input = fread(fid,[nt nx],'single');
fclose(fid);


fid = fopen('receiver_data_p','rb');
primary = fread(fid,[nt nx],'single');
fclose(fid);


for i_iter = 1:nt
  for j_iter = 1:nx
    if ((t(i_iter) < 6.0e-4*x(j_iter) + 0.18) || (t(i_iter) < -6.0e-4*x(j_iter) + 0.18))
      primary(i_iter, j_iter) = 0.0;
    end
  end
end


primary_fk = fft2(primary);


fid = fopen('primary.bin','wb');
fwrite(fid,primary,'single');
fclose(fid);

fid = fopen('primary_fk.bin','wb');
fwrite(fid,abs(primary_fk),'single');
fclose(fid);



exinput = zeros(exnt, exnx);
exinput(:,npad+1:npad+nx) = input;


finput = fft(input);
input_fk = fft2(input);
exfinput = fft(exinput);
exinput_fk = fft2(exinput);


fid = fopen('exdata_with_ghost.bin','w');
fwrite(fid,exinput,'single');
fclose(fid);


fid = fopen('fdata_with_ghost.bin','w');
fwrite(fid,abs(finput),'single');
fclose(fid);


fid = fopen('exfdata_with_ghost.bin','w');
fwrite(fid,abs(exfinput),'single');
fclose(fid);


fid = fopen('data_with_ghost_fk.bin','w');
fwrite(fid,abs(input_fk),'single');
fclose(fid);


fid = fopen('exdata_with_ghost_fk.bin','w');
fwrite(fid,abs(exinput_fk),'single');
fclose(fid);



%% TauP transform
A = zeros(exnp,exnx);
output = zeros(exnt,exnp);
foutput = zeros(exnf,exnp);
tmp_in = zeros(exnx,1);
tmp_out = zeros(exnp,1);

tic

for i_iter = 1:floor(nf/2)
  
  tmp_in = transpose(exfinput(i_iter,:));
    
  for j_iter = 1:exnp
    for k_iter = 1:exnx
      A(j_iter,k_iter) = exp(1i*2.0*pi*exf(i_iter)*exx(k_iter)*exxp(j_iter));   
    end
  end
  
  tmp_out = A*tmp_in;
  
  foutput(i_iter,:) = transpose(tmp_out);
  
end


for i_iter = 2:floor(nf/2)
  foutput(nf-i_iter+2,:) = conj(foutput(i_iter,:)); 
end


output = real(ifft(foutput));

toc



%fdeghost = foutput;
fdeghost = zeros(size(foutput));
tau2_store = zeros(np,1);
tau3_store = zeros(np,1);



%% Filter on every p-trace
for i_iter = 1:np
    
  current_p = i_iter*exdp + exfp;
  
  if (1-(current_p*vwater)^2>0) 
  
    data = output(:,i_iter);
    fdata = foutput(:,i_iter);
    average = mean(abs(data));
    maxvalue = zeros(nt,1);
    testmean = zeros(nt,1);
    check = zeros(nt,1);
    
    for iter = ntcheck_min:ntcheck_max
      tau2 = (iter-1)*dt;
      minustshift = exp(1i*tau2*exomega);
      F2 = 1 - minustshift;
      denominator = 2 -  2*cos(exomega*tau2);
      denominator = denominator + 2.0;
      fP =( fdata.*F2)./denominator;
      for i = nt:-1:floor(nt/2)
        fP(i) = conj(fP(nt - i + 2));
      end

      P2 = real(ifft(fP));
      G2 = data - P2;
 
      maxvalue(iter) = max(P2)/average;
      test = G2 + circshift(P2,[iter-1 0]);
      testmean(iter) = 1./mean(abs(test)/average);

    end

    for i = 1:nt
      check(i) = testmean(i)+maxvalue(i);
    end
  
    [tmp ind] = max(check);

  
    tau2 = (ind -1)*dt;
    tau2_store(i_iter) = tau2;
    tau3 = 2*d*sqrt(1-(current_p*vwater)^2)/vwater;
    tau3_store(i_iter) = tau3;
  
    tau2 = tau3;
  
    minustshift = exp(1i*tau2*exomega);
    F2 = 1 - minustshift;
    denominator = 2 -  2*cos(exomega*tau2);
    denominator = denominator + recover_eps + 20.0*recover_eps*abs(current_p)/pmax;
    ftmp =( fdata.*F2)./denominator;
    for i = nt:-1:floor(nt/2)
      ftmp(i) = conj(ftmp(nt - i + 2));
    end
    fdeghost(:,i_iter) = ftmp;
  
  end
  
end

deghost = real(ifft(fdeghost));



%% Output
fid = fopen('freceiver_data_with_ghost_n1_500_taup.bin','wb');
fwrite(fid,abs(foutput),'single');
fclose(fid);

fid = fopen('receiver_data_with_ghost_n1_500_taup.bin','wb');
fwrite(fid,output,'single');
fclose(fid);


fid = fopen('freceiver_data_with_ghost_n1_500_taup_deghost.bin','wb');
fwrite(fid,abs(fdeghost),'single');
fclose(fid);

fid = fopen('receiver_data_with_ghost_n1_500_taup_deghost.bin','wb');
fwrite(fid,deghost,'single');
fclose(fid);


plot(p,tau2_store,p,real(tau3_store));






