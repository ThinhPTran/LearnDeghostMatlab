clear all
close all
clc


recover_eps=1.0;
nt=1001;
dt=1.3307e-06;
t=(0:dt:(nt-1)*dt)';


fid=fopen('ttrace.bin','rb');
input=fread(fid,[1 nt],'single');
fclose(fid);


figure();
plot(t,input);
legend('input');


df=1.0/(nt*dt);
f=(0:df:(nt-1)*df);
omega=2.0*pi*f;


%% Identify Tau
ntshift=85
tau=ntshift*dt


%% Check tau shift
ntcheck=100;
tcheck=(50*dt:dt:149*dt)';


data = input;
fdata = fft(data);
maxvalue = zeros(ntcheck,1);
testmean = zeros(ntcheck,1);
check = zeros(ntcheck,1);

average = mean(abs(data));

for iter = 1:ntcheck
    tau = tcheck(iter);
    minustshift = exp(1i*tau*omega);
    F2 = 1 - minustshift;
    denominator = 2-2*cos(omega*tau);
    denominator = denominator + 2.0;
    fP =( fdata.*F2)./denominator;
    for i = nt:-1:floor(nt/2)
        fP(i) = conj(fP(nt - i + 2));
    end
    
    P2 = real(ifft(fP));
    G2 = data - P2;
    
    
    maxvalue(iter) = max(P2)/average;
    test = G2 + circshift(P2,[iter-1 0]);
    testmean(iter) = average/mean(abs(test));
    check(iter) = maxvalue(iter)+testmean(iter);
    
end


[tmp ind] = min(check);
ind
tau=tcheck(ind)



figure();
plot(tcheck,check);
legend('check');

  


   
%% Calculate recovered data
data=input;    
fdata=fft(data);
ftmp=zeros(size(fdata));
    

minustshift=exp(1i*tau*omega);
F2=1-minustshift;
denominator=2.0-2.0*cos(omega*tau)+recover_eps;
    
for j_iter = 1:1
    
  minustshift = exp(1i*tau*omega);
  F2 = 1 - minustshift;
  denominator = 2.0 - 2.0*cos(omega*tau)+recover_eps;
  ftmp =( fdata.*F2 + recover_eps*ftmp)./denominator;
  for i = nt:-1:floor(nt/2)
    ftmp(i) = conj(ftmp(nt - i + 2));
  end
    
end
  

output=real(ifft(ftmp));


figure();
plot(t,input,t,output);
legend('input','output');
    







