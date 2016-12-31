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

input=input+2.0*randn(size(input)); 

input1=zeros(size(input)); 


nt=201;
dt = 2e-3;
t=[0:dt:(nt-1)*dt];
df=1./(nt*dt);
f=[0:df:(nt-1)*df];
omega=2.0*pi*f;

hnf=floor(nt/2);


tau=2.0*5.0/1545.8
%tau=2.0*5.0/1500.0
%tau=2.0*5.0/1000.0
eps=0.2


for i_iter = 1:nt
    if ((t(i_iter)>0.19)&&(t(i_iter)<0.198))
       input1(i_iter)=input(i_iter);  
    end
end


figure();
plot(t,input,t,input1);
legend('input','input1');


finput=fft(input);


fabsinput=mag2db(abs(fft(input)));
fabsinput1=mag2db(abs(fft(input1))); 


figure(); 
plot(f,fabsinput,f,fabsinput1); 
legend('fabsinput','fabsinput1'); 




%% Deghost
% First way
nominator1=1-exp(1i*omega*tau);
denominator1=2.0-2.0*cos(omega*tau)+eps;
foutput1=finput.*nominator1./denominator1;

% Second way
output2 = zeros(size(input));

for j_iter = 1:30
    
   foutput2 = fft(output2);

   minustshift = exp(1i*tau*omega);
   nominator2 = 1 - minustshift;
   denominator2 = 2 -  2*cos(omega*tau) + eps;
   finput = fft(input);
   foutput2 =( finput.*nominator2 + eps*foutput2)./denominator2;
   for j = nt:-1:floor(nt/2)
    foutput2(j) = conj(foutput2(nt - j + 2));
   end
 

   output2 = real(ifft(foutput2));
 
end

% Third Way
nominator3=1-exp(1i*omega*tau);
denominator3=2.0-2.0*cos(omega*tau)+eps;
foutput3=finput.*nominator3./denominator3;
foutput3(floor(nt/2)+1:nt)=conj(foutput3(nt+2-[floor(nt/2)+1:nt]));
output3=real(ifft(foutput3)); 
output3=(input+output3)/2.0; 


for i_iter=1:nt 
    if ((abs(output3(i_iter))>0.7*abs(input(i_iter)))&&(abs(output3(i_iter))<1.5*abs(input(i_iter))))
        output3(i_iter)=input(i_iter);
    else
        output3(i_iter)=0; 
    end
end


foutput3=fft(output3); 


% Forth way
nominator4=1-exp(1i*omega*tau);
denominator4=2.0-2.0*cos(omega*tau)+eps;
foutput4=finput.*nominator4./denominator4;
abssmoothfinput=transpose(smooth(abs(finput),30));
angoutput4=angle(foutput4); 
foutput4=abssmoothfinput.*exp(1i*angoutput4); 

% Fifth way
nominator5=1-exp(1i*omega*tau);
denominator5=2.0-2.0*cos(omega*tau)+eps;
foutput5=finput.*nominator5./denominator5;


foutput5(floor(nt/2)+1:nt)=conj(foutput5(nt+2-[floor(nt/2)+1:nt]));
output5=real(ifft(foutput5));

input2=(input+output5)/2.0;

for i_iter=1:nt 
    if ((abs(input2(i_iter))>0.7*abs(input(i_iter)))&&(abs(input2(i_iter))<1.5*abs(input(i_iter))))
        input2(i_iter)=input(i_iter);
    else
        input2(i_iter)=0; 
    end
end

finput2=fft(input2);


figure();
plot(t,input,t,input2,t,output5)
legend('input','input2','output5')


abssmoothfinput=transpose(smooth(abs(finput2)));
angoutput5=angle(foutput5); 
foutput5=abssmoothfinput.*exp(1i*angoutput5); 




ffilter1=abs(nominator1./denominator1);
ffilter2=abs(nominator2./denominator2);



% for i=nt:-1:floor(nt/2)
%    foutput(i)=conj(foutput(nt-i+2)); 
% end


foutput1(floor(nt/2)+1:nt)=conj(foutput1(nt+2-[floor(nt/2)+1:nt]));
foutput2(floor(nt/2)+1:nt)=conj(foutput2(nt+2-[floor(nt/2)+1:nt]));
foutput3(floor(nt/2)+1:nt)=conj(foutput3(nt+2-[floor(nt/2)+1:nt]));
foutput4(floor(nt/2)+1:nt)=conj(foutput4(nt+2-[floor(nt/2)+1:nt]));
foutput5(floor(nt/2)+1:nt)=conj(foutput5(nt+2-[floor(nt/2)+1:nt]));


output1=real(ifft(foutput1));
output2=real(ifft(foutput2));
output3=real(ifft(foutput3));
output4=real(ifft(foutput4)); 
output5=real(ifft(foutput5)); 



%% Mute both sides for a clean output
% [max_out1,nmax1]= max(output1);
% [max_out2,nmax2]= max(output2);
% 
% output_cut1=[output1(nmax1-25:nmax1+25)];  
% output_cut1([1:10 end-9:end])=0;
% output_sym1 =[output1(nmax1-25:nmax1) output1(nmax1-1:-1:nmax1-25)];
% output_sym1([1:10 end-9:end])=0;


figure();
plot(t,input,t,output1,t,output2,t,output3,t,output4,t,output5);
% plot(t,input,t,output3); 
legend('input','deghost1','deghost2','deghost3','deghost4','deghost5');


% figure();
% plot(t,input,t,output3)
% legend('input','deghost3'); 


%% Check the spectrum 
finput=mag2db(abs(fft(input)));
foutput1=mag2db(abs(fft(output1)));
foutput2=mag2db(abs(fft(output2)));
foutput3=mag2db(abs(fft(output3))); 
foutput4=mag2db(abs(fft(output4))); 
foutput5=mag2db(abs(fft(output5))); 


figure();
plot(f(1:hnf),finput(1:hnf),f(1:hnf),foutput1(1:hnf),f(1:hnf),foutput2(1:hnf),f(1:hnf),foutput3(1:hnf),f(1:hnf),foutput4(1:hnf),f(1:hnf),foutput5(1:hnf));
legend('finput','foutput1','foutput2','foutput3','foutput4','foutput5');


%% Check the filter spectrum
figure();
plot(f(1:hnf),ffilter1(1:hnf),f(1:hnf),ffilter2(1:hnf));
legend('filter 1','filter 2');



% figure;
% plot(0:dt:(length(output_sym)-1)*dt,output_sym,0:dt:(length(output_sym)-1)*dt,output_cut);
% 
% 
% figure,
% plot([1:length(t)]/max(t),abs(real(fft(input))),'b',[1:length(t)]/max(t),abs(real(fft(output))),'r');
% legend('input','output');


% save('farfield_wavelet_sym.txt','output_sym1','-ascii');
% save('farfield_wavelet.txt','output_cut1','-ascii');


%%% Finer grid source 
% finput=fft(output);
% finedt=dt/20.0;
% finet= 0:finedt:max(t);
% finent = length(finet);
% 
% finefinput=zeros(size(finet));
% finefinput(1:floor(nt/2)-1)=finput(1:floor(nt/2)-1);
% 
% for i_iter=finent:-1:floor(finent/2)
%   finefinput(i_iter)=conj(finefinput(finent-i_iter+2)); 
% end
% 
% fineinput=real(ifft(finefinput));


% figure();
% plot(t,output);
% legend('output');


% figure();
% plot(finet,fineinput);
% legend('fineinput');


% fid=fopen('farfield_sig','w');
% fwrite(fid,input,'single');
% fclose(fid);
% 
% 
% figure();
% plot(t,input);
% legend('input');














