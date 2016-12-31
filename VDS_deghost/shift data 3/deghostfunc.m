function [output zout] = deghostfunc(input, zin, vwater, eps, nt, dt)

%% Parameter
% Basic params
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


zmin_check = max(1,zin-2);
zmax_check = zin+2;


dt_check = dt/4.0;
tmin_check = floor((2.0*zmin_check/vwater)/dt_check)*dt_check;
tmax_check = floor((2.0*zmax_check/vwater)/dt_check)*dt_check;

t_check = (tmin_check:dt_check:tmax_check)';
nt_check = length(t_check);


tau = 2.0*zin/vwater;


% Clean-up params
fmin=0.01;
fmax=(nt-2)*df;
res=4.0;
nfb=15;



%% Calculate norms
finput = fft(input); 


normtestmean = zeros(nt_check,1);
normmaxvalue = zeros(nt_check,1);
normderiv = zeros(nt_check,1);
normderiv2 = zeros(nt_check,1);
normcheck = zeros(nt_check,1);
normtest = zeros(nt_check,1);


average = mean(abs(input));


for iter = 1:nt_check
    
 newP = zeros(size(input));
 P = zeros(size(input));

 ctau = t_check(iter);
 
 for j_iter = 1:30
    
   P = newP; 
   fP = fft(P);

   minustshift = exp(1i*ctau*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*ctau) + eps;
   finput = fft(input);
   fP =( finput.*F2 + eps*fP)./denominator;
   for j = nt:-1:floor(nt/2)
    fP(j) = conj(fP(nt - j + 2));
   end
 

   newP = real(ifft(fP));
   diff = abs(newP - P);
   
   f = mean(abs(diff))/mean(abs(newP));
   
   if (f < 0.01) 
      break; 
   end
 
 end
 
 j_iter;
 
 G = input - P;
 

 % Calculate normmax value
 normmaxvalue(iter) = max(P)/average;
 
 % Calculate normtestmean value
 fP = fft(P);
 shiftfP = fP.*exp(-1i*ctau*omega);
 for j = nt:-1:floor(nt/2)
    shiftfP(j) = conj(shiftfP(nt - j + 2));
 end
 shiftP = real(ifft(shiftfP));
 test = shiftP+G;
 normtestmean(iter) = average/mean(abs(test));
 
 % Calculate normcheck
 normcheck(iter) = 1.0/mean(abs(P));
 
 
%  figure();
%  plot(t,P,t,shiftP,t,G,t,withghost,t,test);
%  legend('P','shiftP','G','withghost','test');
 
 
end


% Normalize
normmaxvalue=normmaxvalue/(max(abs(normmaxvalue)));
normtestmean=normtestmean/(max(abs(normtestmean)));
normcheck=normcheck/(max(abs(normcheck)));

% Calculate normderiv
normderiv=abs(circshift(normmaxvalue,[1 0]) - circshift(normmaxvalue,[-1 0]));
normderiv2=abs(circshift(normmaxvalue,[1 0]) - 2.0*normmaxvalue + circshift(normmaxvalue,[-1 0]));

normderiv=1.0./normderiv;
normderiv2=1.0./normderiv2;

normderiv=normderiv/(max(abs(normderiv)));
normderiv2=normderiv2/(max(abs(normderiv2)));


% Calculate norm for test
% if (z<=6)
% normtest=(normmaxvalue+2.0*normderiv2+2.0*normtestmean+normcheck)/7.0;
% else 
normtest=(normmaxvalue+normtestmean+normcheck)/3.0;    
% end


% figure()
% plot(t_check,normmaxvalue,t_check,normtestmean,t_check,normtest,t_check,normcheck);
% legend('normmaxvalue','normtestmean','normtest','normcheck');


[tmp ind] = max(normtest);



%% Recover
 est_tau = (ind-1)*dt_check+tmin_check;
 est_tau = tau;
 
 
 P = zeros(size(input));
 
% Method 1
%   minustshift = exp(1i*est_tau*omega);
%   F2 = 1 - minustshift;
%   denominator = 2 -  2*cos(omega*est_tau) + eps;
%   finput = fft(input);
%   fP =( finput.*F2)./denominator;
%   for j = nt:-1:floor(nt/2)
%     fP(j) = conj(fP(nt - j + 2));
%   end
%   
%   P = real(ifft(fP)); 
 
 
% Method 2
 for j_iter = 1:30
    
   fP = fft(P);

   minustshift = exp(1i*est_tau*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*est_tau) + eps;
   finput = fft(input);
   fP =( finput.*F2 + eps*fP)./denominator;
   for j = nt:-1:floor(nt/2)
    fP(j) = conj(fP(nt - j + 2));
   end
 

   P = real(ifft(fP));
 
 end
 
 
 
 %% Stabilize the spect
%  fP = fft(P); 
%  meanfP = mean(abs(fP));
%  
%  
%  for iter = 1:nt
%    if (abs(fP(iter))>2.0*meanfP)
%      tmp = angle(fP(iter)); 
%      fP(iter) = 0.7*meanfP*exp(1i*tmp); 
%    end
%  end
%  
%  
%  P = real(ifft(fP)); 
 
 
 
 %% Clean up process
 P = cleanup(input, P , nt, dt, nfb, fmin, fmax, res);
 
 
 
 %% Output
 output = P; 
 
 zout = est_tau*vwater/2.0;
 







