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
hnf=floor(nt/2);

dt = 2e-3;
t=[0:dt:(nt-1)*dt];
df=1./(nt*dt);
f=[0:hnf -hnf:-1]*df;
omega=2.0*pi*f;


tau=2.0*5.0/1545.8
% tau=0.0007
% tau = 0; 
%tau=2.0*5.0/1500.0
%tau=2.0*5.0/1000.0
eps=0.1

% nominator=1-exp(-1i*omega*0.008);
% 
% finput=fft(input); 
% finput=finput.*nominator;
% input=real(ifft(finput)); 
% 
% input=input+2.0*randn(size(input)); 


figure();
plot(t,input);
title('input'); 


finput=fft(input);


fabsinput=mag2db(abs(fft(input)));


figure(); 
plot(f(1:hnf),fabsinput(1:hnf)); 
legend('fabsinput'); 




%% Deghost
% First way
nominator1=1-exp(1i*omega*tau);
denominator1=2.0-2.0*cos(omega*tau)+eps;
ffilter1=nominator1./denominator1;
filter1=real(ifft(ffilter1));  
absffilter1=abs(ffilter1); 
foutput1=finput.*ffilter1;


% Second way
ffilter2=1.0./(1-exp(-1i*omega*tau)+eps*(1+exp(-1i*omega*tau)));
absffilter2=abs(ffilter2); 
foutput2=finput.*ffilter2; 


% Third way
ffilter3=(1-exp(-1i*omega*tau));
absffilter3=abs(ffilter3); 
angffilter3=angle(ffilter3); 
index=find(absffilter3<=eps); 
absffilter3(index)=eps;
ffilter3=absffilter3.*exp(1i*angffilter3); 
ffilter3=1.0./ffilter3; 
absffilter3=abs(ffilter3); 
foutput3=finput.*ffilter3; 


% Forth way
nominator4=1-exp(1i*omega*tau);
denominator4=2.0-2.0*cos(omega*tau)+eps;
ffilter4=nominator4./denominator4;
absffilter4=abs(ffilter4); 
angffilter4=smooth(angle(ffilter4),30)'; 
ffilter4=absffilter4.*exp(1i*angffilter4); 
foutput4=finput.*ffilter4;


% Fifth way
nominator5=1-exp(1i*omega*tau); 
denominator5=2.0-2.0*cos(omega*tau)+eps;
ffilter5=nominator5./denominator5; 
filter5=real(ifft(ffilter5)); 
fillngth=0.2; 
index=find(t>fillngth); 
filter5(index)=0; 
for i_iter = 2:hnf
   filter5(nt-i_iter+2)=-filter5(i_iter);  
end
ffilter5=fft(filter5); 
foutput5=finput.*ffilter5; 


% Sixth way
nominator6=1-exp(1i*omega*tau); 
denominator6=2.0-2.0*cos(omega*tau)+eps; 
ffilter6=nominator6./denominator6;
%ffilter6=ffilter3; 
filter6=real(ifft(ffilter6)); 
c=fillngth/2.00; 
tapper=0.3*fillngth/2.0;
gauw=1.2*exp(-(t.*t)/(2.0*c*c)); 
linear=(1.0-(t-fillngth-tapper)/(tapper));
index=find(t<fillngth-tapper); 
gauw(index)=1.0; 
linear(index)=1.0; 
index=find(t>fillngth); 
linear(index)=0.0;

nfill=floor(fillngth/(2.0*dt)); 
ntap=floor(tapper/(2.0*dt)); 

n1=(ntap:2*ntap);
truncfil1=zeros(1,nt);

truncfil1(1:nfill-ntap)=1.0;
truncfil1(nfill-ntap:nfill)=0.5-0.5*cos(2*pi*n1/(2*ntap)); 
truncfil1(nfill+1:nt)=0.0; 


filter6=filter6.*linear; 
filter7=filter6.*truncfil1;
filter6(hnf:nt)=0; 

for i_iter = 2:hnf
   filter6(nt-i_iter+2)=-filter6(i_iter);  
   filter7(nt-i_iter+2)=-filter7(i_iter); 
end
ffilter6=fft(filter6); 
ffilter7=fft(filter7); 
foutput6=finput.*ffilter6;
foutput7=finput.*ffilter7; 


figure(); 
plot(t,gauw,t,linear); 
legend('gauw','linear'); 


output1=real(ifft(foutput1));
output2=real(ifft(foutput2)); 
output3=real(ifft(foutput3)); 
output4=real(ifft(foutput4)); 
output5=real(ifft(foutput5));
output6=real(ifft(foutput6)); 
output7=real(ifft(foutput7)); 



figure();
plot(t,input,t,output1,t,output2,t,output3); 
legend('input','deghost1','deghost2','deghost3');

figure();
plot(t,input,t,output1,t,output4); 
legend('input','deghost1','deghost4');


figure(); 
plot(t,filter1,t,filter5,t,filter6); 
legend('filter1','filter5','filter6'); 

figure();
plot(t,input,t,output1); 
legend('input','deghost1');

figure(); 
plot(t,input,t,output1,t,output5,t,output6); 
legend('input','deghost1','deghost5','deghost6'); 


% figure();
% plot(t,input,t,output1,t,output4); 
% legend('input','deghost1','deghost4');



%% Check the spectrum 
finput=mag2db(abs(fft(input)));
foutput1=mag2db(abs(fft(output1)));
foutput2=mag2db(abs(fft(output2))); 
foutput3=mag2db(abs(fft(output3))); 
foutput4=mag2db(abs(fft(output4))); 
foutput5=mag2db(abs(fft(output5))); 
foutput6=mag2db(abs(fft(output6))); 
foutput7=mag2db(abs(fft(output7))); 
filter1=mag2db(abs(absffilter1));
filter2=mag2db(abs(absffilter2));
filter3=mag2db(abs(absffilter3));



figure();
plot(f(1:hnf),finput(1:hnf),f(1:hnf),foutput1(1:hnf),f(1:hnf),foutput2(1:hnf),f(1:hnf),foutput3(1:hnf));
legend('finput','foutput1','foutput2','foutput3');

figure();
plot(f(1:hnf),finput(1:hnf));
legend('finput');
ylim([-50 50]); 

figure();
plot(f(1:hnf),foutput1(1:hnf));
legend('foutput1');
ylim([-50 50]); 

figure();
plot(f(1:hnf),foutput5(1:hnf));
legend('foutput5');
ylim([-50 50]); 

figure();
plot(f(1:hnf),finput(1:hnf),f(1:hnf),foutput1(1:hnf));
legend('finput','foutput1');
ylim([-50 50]); 


%% Check the filter spectrum
% figure();
% plot(f(1:hnf),filter1(1:hnf),f(1:hnf),filter2(1:hnf),f(1:hnf),filter3(1:hnf));
% legend('filter 1','filter 2','filter 3');


%% Check phase of filter
% figure()
% plot(f(1:hnf),angle(ffilter1(1:hnf)));
% legend('phase f1'); 
% 
% figure()
% plot(f(1:hnf),angle(ffilter2(1:hnf)));
% legend('phase f2'); 
% 
% figure()
% plot(f(1:hnf),angle(ffilter3(1:hnf)));
% legend('phase f3'); 
% 
% figure()
% plot(f(1:hnf),angle(ffilter4(1:hnf)));
% legend('phase f4'); 







