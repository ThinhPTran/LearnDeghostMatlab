iter_cg = 10000;
epsilon = 0.001;

nt = 2001;
n2 = 1001;

t = (0:0.004:8)';
f = (0:0.125:250)';
omega = 2*pi*f;
z = 80.0;
tau = 2.0*z/1500;
tshift = exp(-1i*tau*omega);
F = 1 .- tshift;


fid = fopen('primary.bin','r');
synthetic = fread(fid,[nt n2],'float');
fclose(fid);

%%Get primary
primary = synthetic(1:nt,floor(n2/2));
primary += 1./cosh(100*(t-2)) + 1./cosh(100*(t-4));
fprimary = fft(primary);
figure(1)
plot(t,primary);
%

%%Synthetize data with ghost
fwithghost = fprimary.*F;
fwithghost(1) = 1.0;
for i = 2001:-1:floor(2001/2)
  fwithghost(i) = conj(fwithghost(2001 - i + 2));
endfor
tmp= ifft(fwithghost);
noise = randn(2001,1)/10.0;
withghost = real(tmp) + noise;
fwithghost = fft(withghost);
figure(2)
plot(t,real(withghost));
%


%%Find the Primary (Deghosting) 

%%Initial value of tau
tau = 0.1066666;
t1shift = exp(-1i*tau*omega);
F1 = 1 .- t1shift;
for i = 2001:-1:floor(2001/2)
F1(i) = F1(2001 - i + 2);	
endfor
%%%%%%%%%%%%%%%%%%%%%


%%Calculate Primary 
for i = 1:10
  	
%if ( i == 1 )
  W = abs(fwithghost);
%else 
  %W = abs(absfP);
%endif

X = cal_fprimary(F1,W,fwithghost, iter_cg, epsilon);

for j = 2001:-1:floor(2001/2)
  X(j) = conj(X(2001 - j + 2));
endfor

X(1) = 1.0;
absftrace = abs(fwithghost);
absfP = abs(X);
ratio = absftrace./absfP;

adjustfP = X.*ratio;

tmp = ifft(adjustfP);
P = real(tmp);

figure(3);
plot(t,P,'color','red');

G = withghost - P;
figure(4);
plot(t,G,'color','red');
fG = fft(G);

%Find exp(1i*tau*omega);
X = cal_exptau(fG,-adjustfP,iter_cg,epsilon);
t1shift = conj(X);
F1 = 1 .- t1shift;
for j = 2001:-1:floor(2001/2)
F1(j) = F1(2001 - j + 2);	
endfor

endfor
% End calculate Primary


figure(5);
plot(t,primary,'color','blue',t,withghost,'color','red',t,P,'color','green',t,G,'color','yellow');
%










