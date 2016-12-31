clear all
close all
clc

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
d = 20;
v = 1500;
zmin_check = 5;
zmax_check = 200;

taumin_check = 2.0*zmin_check/v;
taumax_check = 2.0*zmax_check/v;

ntcheck_min = floor(taumin_check/dt)+1;
ntcheck_max = floor(taumax_check/dt)+1;




%% input
fid = fopen('data_with_ghost','r');
input = fread(fid,[nt nx],'single');
fclose(fid);


finput = fft(input);

fid = fopen('fdata_with_ghost','w');
fwrite(fid,abs(finput),'single');
fclose(fid);



%% FK transform
input_fk = fft2(input);

fdeghost = zeros(size(input_fk));


%% Filter on every k-trace
for i_iter = 1:floor(nkx/2)
  
  fdata = input_fk(:,i_iter);
  data = ifft(fdata);
  average = mean(abs(data));
  maxvalue = zeros(nt,1);
  testmean = zeros(nt,1);
  check = zeros(nt,1);
    
  for iter = ntcheck_min:ntcheck_max
    tau2 = (iter-1)*dt;
    minustshift = exp(1i*tau2*omega);
    F2 = 1 - minustshift;
    denominator = 2 -  2*cos(omega*tau2);
    denominator = denominator + 2.0;
    fP =( fdata.*F2)./denominator;

    P2 = ifft(fP);
    G2 = data - P2;
 
    maxvalue(iter) = max(abs(P2))/average;
    test = G2 + circshift(P2,[iter-1 0]);
    testmean(iter) = 1./mean(abs(test)/average);

  end

  for i = 1:nt
    check(i) = testmean(i)+maxvalue(i);
  end
  
  [tmp ind] = max(check);

  current_kx = i_iter*dkx
  tau2 = (ind -1)*dt 
  
  minustshift = exp(1i*tau2*omega);
  F2 = 1 - minustshift;
  denominator = 2 -  2*cos(omega*tau2);
  denominator = denominator + 0.1;
  ftmp =( fdata.*F2)./denominator;
  fdeghost(:,i_iter) = ftmp;
end


deghost = real(ifft2(fdeghost));



%% Output
fid = fopen('deghost.bin','w');
fwrite(fid,deghost,'single');
fclose(fid);



