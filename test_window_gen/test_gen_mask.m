clear all
close all
clc


nt=1001;
dt=0.004;
t=(0:dt:(nt-1)*dt);
safegap=0.02;
ns=floor(safegap/(2.0*dt)); 


data=zeros(1,nt); 

data(100)=1;
data(105)=1;
data(200)=1;
data(500)=1; 
data(800)=1; 


squarepulse=zeros(nt,1); 
squarepulse(1:6)=1;
squarepulse(997:nt)=1; 
gausspulse=1.2*exp(-(t.*t)/(2.0*(safegap/2.3)*(safegap/2.3)));
cosine=cos(pi*t/(nt*dt)-pi/2.0); 
cossqr=cosine.*cosine; 


gausspulse(1:ns+1)=1.0; 


for i_iter = 2:floor(nt/2)
    gausspulse(nt-i_iter+2)=gausspulse(i_iter);
end



figure(); 
plot(t,squarepulse); 
title('squarepulse'); 


figure();
plot(t,gausspulse); 
title('gausspulse'); 


figure();
plot(t,cosine); 
title('cosine'); 


figure();
plot(t,cossqr); 
title('cossqr'); 


figure(); 
plot(t,data); 
title('data'); 


fsquarepulse=fft(squarepulse); 
fdata=fft(data);

fgausspulse=fft(gausspulse); 


foutput=fgausspulse.*fdata; 

output=real(ifft(foutput));


figure(); 
plot(t,output,t,data); 
legend('output','data'); 





