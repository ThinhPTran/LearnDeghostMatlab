clear all
close all
clc


nt=201;
nh=floor(nt/2); 
n=nh:nt-1;
dt=0.004; 

df=1/(nt*dt); 
f=[0:nt-1]*df; 


size(n)

alpha=0.5;
beta=0.5;


w=alpha-beta*cos(2*pi*n/(nt-1)); 


figure(); 
plot(n,w); 
legend('w'); 


nfill=25; 
ntap=10; 
n1=(ntap:2*ntap)';
n2=0:ntap; 

truncfil1=zeros(nt,1);
truncfil2=zeros(nt,2); 

truncfil1(1:nfill-ntap)=1.0;
truncfil2(1:nfill-ntap)=1.0; 
truncfil1(nfill-ntap:nfill)=alpha-beta*cos(2*pi*n1/(2*ntap)); 
truncfil2(nfill-ntap:nfill)=1.0-n2/ntap; 
truncfil1(nfill+1:nt)=0.0; 
truncfil2(nfill+1:nt)=0.0;

for i_iter=2:nh
   truncfil1(nt-i_iter+2) = truncfil1(i_iter); 
end


figure(); 
plot(1:nt,truncfil1,1:nt,truncfil2); 
legend('truncfil1','truncfil2'); 


ffil1=fft(truncfil1);
ffil2=fft(truncfil2); 


figure(); 
plot(f(1:nh),mag2db(abs(ffil1(1:nh)))); 
legend('ffil1'); 











