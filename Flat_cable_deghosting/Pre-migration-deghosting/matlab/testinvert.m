nt = 2001;
n2 = 1001;


z = 50.0;
tau = 2.0*z/1500;
tshift = exp(-1i*tau*omega);
F = 1 - tshift;
F(1) = 1.0;

fid = fopen('primary.bin','r');
synthetic = fread(fid,[nt n2],'float');
fclose(fid);

%%Get primary
tmp1 = synthetic(1:nt,floor(n2/2));
tmp2 = 1./cosh(100*(t-2)) + 1./cosh(100*(t-4));
primary_sig = tmp1 + tmp2;
fprimary_sig = fft(primary_sig);
figure(1)
plot(t,primary_sig);
%%%%%


%%Synthetize data with ghost
fwithghost = fprimary_sig.*F;
fwithghost(1) = 1.0;
for i = 2001:-1:floor(2001/2)
  fwithghost(i) = conj(fwithghost(2001 - i + 2));
end
tmp= ifft(fwithghost);
noise = randn(2001,1)/10.0;
withghost = real(tmp);% + noise;
fwithghost = fft(withghost);
figure(2)
plot(t,real(withghost));
%%

%% Calculate Primary (Method I)
z1 = 50;
tau1 = 2.0*z1/1500.0;
tshift1 = exp(-1i*tau1*omega);
F1 = 1 - tshift1;
F1(1) = 1.0;
fP = fwithghost./F1;
for i = 2001:-1:floor(2001/2)
  fP(i) = conj(fP(2001 - i + 2));
end

%absftrace = abs(fwithghost);
%absfP = abs(fP);
%ratio = absftrace./absfP;
%adjustfP = fP.*ratio;

P1 = ifft(fP);
figure(3)
plot(t,real(P1));
%%

%% Calculate Primary (Method II)
z2 = 50;
tau2 = 2.0*z2/1500.0;
minustshift = exp(1i*tau2*omega);
F1 = 1 - minustshift;
denominator = 2 -  2*cos(omega*tau2) + 1e-15;
fP =( fwithghost.*F1)./denominator;
for i = 2001:-1:floor(2001/2)
  fP(i) = conj(fP(2001 - i + 2));
end

P2 = ifft(fP);
figure(4)
plot(t,real(P2));
%% 

%% Calculate Primary (Method III)


%%
