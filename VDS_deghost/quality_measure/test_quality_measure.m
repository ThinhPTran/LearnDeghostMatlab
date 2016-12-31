clc
clear all;
close all;

input=[0.000     0.012     0.032     0.044     0.020    -0.026    -0.068 ...
-0.077    -0.040     0.010     0.057     0.076     0.045    -0.004 ...
-0.051    -0.075    -0.052    -0.005     0.043     0.075     0.059 ...
 0.012    -0.037    -0.072    -0.065    -0.023     0.026     0.068 ...
 0.071     0.032    -0.018    -0.062    -0.073    -0.039     0.010 ...
 0.057     0.075     0.046    -0.003    -0.050    -0.074    -0.052 ...
-0.004     0.045     0.077     0.061     0.014    -0.036    -0.073 ...
-0.066    -0.023     0.027     0.070     0.072     0.032    -0.016 ...
-0.061    -0.073    -0.039     0.009     0.055     0.073     0.044 ...
-0.004    -0.051    -0.076    -0.056    -0.009     0.053     0.081 ...
 0.069     0.017    -0.055    -0.048    -0.104     0.013     0.002 ...
 0.057     0.129    -0.096     0.165    -0.291     0.153    -0.198 ...
 0.033     0.250    -0.383     0.783    -0.972     1.035    -1.088 ...
 0.623    -0.024    -0.995     2.703    -5.253    30.984    50.263 ...
22.011   -16.056   -59.076   -33.010    -7.162    -8.071    -1.255 ...
-6.592    -0.126     0.063    -0.782     1.831    -0.028     1.655 ...
 0.826     1.223     1.453     1.137     1.528     1.381     1.578 ...
 1.615     1.601     1.594     1.538     1.519     1.548     1.593 ...
 1.670     1.716     1.734     1.722     1.689     1.660     1.629 ...
 1.626     1.609     1.569     1.508     1.416     1.333     1.234 ...
 1.073     0.810     0.460     0.165     0.053     0.125     0.297 ...
 0.470     0.584     0.615     0.531     0.332     0.050    -0.331 ...
-0.840    -1.405    -1.866    -2.047    -1.886    -1.504    -1.050 ...
-0.671    -0.510    -0.641    -1.036    -1.470    -1.664    -1.530 ...
-1.196    -0.854    -0.609    -0.460    -0.387    -0.396    -0.482 ...
-0.614    -0.744    -0.807    -0.795    -0.749    -0.691    -0.640 ...
-0.560    -0.414    -0.231    -0.030     0.139     0.244     0.310 ...
 0.352     0.400     0.476     0.560     0.629     0.668     0.651 ...
 0.598     0.539     0.492     0.481     0.490];

%input(floor(3.0*end/4.0):end)=0;

nt=201;
nf=nt; 
nh=floor(nt/2); 
dt = 2e-3;
t=[0:dt:(nt-1)*dt];
df=1./(nt*dt);
f=[0:nh -nh:-1]*df;
omega=2.0*pi*f;

hnf=floor(nt/2);

vwater=1545.8;
tau=2.0*5.0/1545.8
% tau=2.0*5.0/1000.0
eps=0.0002


input=input+0.0*randn(1,nt); 


figure();
plot(t,input);
legend('input');


finput=fft(input);


%% Deghost
% First way
nominator=1-exp(1i*omega*tau);
denominator=2.0-2.0*cos(omega*tau)+eps;
% foutput=finput.*nominator./denominator;


nominator1=1-exp(1i*omega*tau); 
nominator2=1+exp(1i*omega*tau); 

foutput1=finput.*nominator1;
foutput2=finput.*nominator2; 


output1=real(ifft(foutput1));
output2=real(ifft(foutput2)); 

for i_iter = 1:nt
  if ((output1(i_iter)*output2(i_iter))<0)
    output3(i_iter) = output1(i_iter) + output2(i_iter);
  else
    output3(i_iter) = output1(i_iter) - output2(i_iter); 
  end
end

output3=output3/2.0; 


figure();
plot(t,input,t,output1,t,output2,t,output3);
legend('input','output1','output2','output3');



absfinput=abs(fft(input)); 
absfoutput1=abs(fft(output1)); 
absfoutput2=abs(fft(output2)); 
absfoutput3=abs(fft(output3)); 
dbabsfinput=mag2db(absfinput); 
dbabsfoutput1=mag2db(absfoutput1); 
dbabsfoutput2=mag2db(absfoutput2); 
dbabsfoutput3=mag2db(absfoutput3); 


figure();
plot(1:nf,absfinput,1:nf,absfoutput1,1:nf,absfoutput2,1:nf,absfoutput3); 
legend('absfinput','absfoutput1','absfoutput2','absfoutput3'); 

figure(); 
plot(1:nf,dbabsfinput,1:nf,dbabsfoutput1,1:nf,dbabsfoutput2,1:nf,dbabsfoutput3); 
legend('dbabsfinput','dbabsfoutput1','dbabsfoutput2','dbabsfoutput3'); 



%% Check quality measure
dtau=0.0005; 
ntest=21; 

qualmes=zeros(ntest,1); 
qualmes1=zeros(ntest,1); 
sumabs=zeros(ntest,1); 
derive=zeros(ntest,1); 



for i_iter=1:ntest
    
   ctau=tau+(i_iter-floor(ntest/2)-1)*dtau;
   
   nominator=1-exp(1i*omega*ctau);
   denominator=2.0-2.0*cos(omega*ctau)+eps;
   cfoutput=finput.*nominator./denominator;
   
   cabsfoutput=abs(cfoutput); 
   cangfoutput=angle(cfoutput); 
   
   cabsfoutput=transpose(smooth(abs(finput),30)); 
   cabsfoutput1=abs(finput); 
   
   cfoutput=cabsfoutput.*exp(1i*cangfoutput); 
   cfoutput1=cabsfoutput1.*exp(1i*cangfoutput); 

   coutput=real(ifft(cfoutput)); 
   coutput1=real(ifft(cfoutput1)); 
   
   qualmes(i_iter)=max(coutput); 
   qualmes1(i_iter)=max(coutput1); 
   sumabs(i_iter)=sum(abs(coutput).*abs(coutput)); 
   derive(i_iter)=1./sum(abs(diff(coutput))); 

    
end


% Normalize 
qualmes=qualmes/(max(qualmes));
qualmes1=qualmes1/(max(qualmes1));
derive=derive/(max(derive));



figure(); 
plot(1:ntest,qualmes,1:ntest,qualmes1); 
legend('qualmes','qualmes1'); 

figure(); 
plot(1:ntest,derive); 
title('derivative'); 

figure();
plot(1:ntest,sumabs); 
title('sumabs'); 




[maxval,maxidx]=max(qualmes);

maxidx

ctau=tau+(maxidx-floor(ntest/2)-1)*dtau
zest=ctau*vwater/2.0



%% Check the spectrum 
% finput=mag2db(abs(fft(input)));
% foutput=mag2db(abs(fft(output)));
% 
% 
% figure();
% plot(f(1:hnf),finput(1:hnf),f(1:hnf),foutput(1:hnf));
% legend('finput','foutput');

















