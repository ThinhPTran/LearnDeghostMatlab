function output = deghostfunc(input, zrecin, zsrcin, vwater, eps, nt, dt)

%% Parameter
% Basic params
df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;



rectau = 2.0*zrecin/vwater;
srctau = 2.0*zsrcin/vwater;


% Clean-up params
fmin=0.01;
fmax=(nt-2)*df;
res=4.0;
nfb=35;



%% Recover 
P = zeros(size(input));
shiftinput = zeros(size(input)); 
finput = fft(input); 
 

% Calculate Pspec
Pspec = zeros(size(input)); 
Pspec1 = zeros(size(input)); 


% Method 2
% De recevier ghost
 for j_iter = 1:30
    
   fPspec1 = fft(Pspec1);

   minustshift = exp(1i*rectau*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*rectau) + eps;
   finput = fft(input);
   fPspec1 =( finput.*F2 + eps*fPspec1)./denominator;
   for j = nt:-1:floor(nt/2)
    fPspec1(j) = conj(fPspec1(nt - j + 2));
   end
 

   Pspec1 = real(ifft(fPspec1));
 
 end
 
 
 % De src ghost
 for j_iter = 1:30
    
   fPspec = fft(Pspec);

   minustshift = exp(1i*srctau*omega);
   F2 = 1 - minustshift;
   denominator = 2 -  2*cos(omega*srctau) + eps;
   fPspec1 = fft(Pspec1);
   fPspec =( fPspec1.*F2 + eps*fPspec)./denominator;
   for j = nt:-1:floor(nt/2)
    fPspec(j) = conj(fPspec(nt - j + 2));
   end
 

   Pspec = real(ifft(fPspec));
 
 end
 
 
Pspec = cleanup(input, Pspec, nt, dt, nfb, fmin, fmax, res); 



% Calculate P1
fshiftinput = finput.*exp(1i*rectau*omega); 
for j = nt:-1:floor(nt/2)
   fshiftinput(j) = conj(fshiftinput(nt - j + 2));  
end
shiftinput = -real(ifft(fshiftinput)); 

% P1orig = shiftinput; 
P1 = cleanup(input, shiftinput, nt, dt, nfb, fmin, fmax, res);


% Calculate P2
fshiftinput = finput.*exp(1i*srctau*omega); 
for j = nt:-1:floor(nt/2)
   fshiftinput(j) = conj(fshiftinput(nt - j + 2));  
end
shiftinput = -real(ifft(fshiftinput)); 

% P2orig = shiftinput; 
P2 = cleanup(input, shiftinput, nt, dt, nfb, fmin, fmax, res); 


% Calculate P3
fshiftinput = finput.*exp(1i*(rectau+srctau)*omega); 
for j = nt:-1:floor(nt/2)
   fshiftinput(j) = conj(fshiftinput(nt - j + 2));  
end
shiftinput = real(ifft(fshiftinput));

% P3orig = shiftinput; 
P3 = cleanup(input, shiftinput, nt, dt, nfb, fmin, fmax, res); 



% Plot to see
% Porig = (input + P1orig + P2orig + P3orig)/2.0; 


% figure(); 
% plot(t,input,t,P1orig,t,P2orig,t,P3orig,t,Porig); 
% legend('input','P1orig','P2orig','P3orig','Porig'); 
% title('before clean'); 

figure(); 
plot(t,input,t,P1,t,P2,t,P3,t,Pspec); 
legend('input','P1','P2','P3','Pspec'); 
title('after clean'); 


% Get the coherence
[input P1 P2 P3] = getcoherence(input, P1, P2, P3, nt, dt, nfb, fmin, fmax, res);
% Porigclean = cleanup(input, Porig, nt, dt, nfb, fmin, fmax, res); 

 
%P = (P1+P2);
P = (input+P1+P2+P3)/2.0; 


% Scale back 
maxinput = max(input); 
maxP = max(P); 

P = P*maxinput/maxP; 


% Plot to see
figure(); 
plot(t,input,t,P1,t,P2,t,P3,t,Pspec,t,P); 
legend('input','P1','P2','P3','Pspec','P'); 
title('Before improve'); 


% Improve effect
Pout = improveeffect(input, P, nt);


% P = (P + Pout)/2.0; 
P = Pout;  


% Plot to see
figure(); 
plot(t,input,t,P1,t,P2,t,P3,t,Pspec,t,P); 
legend('input','P1','P2','P3','Pspec','P'); 
title('After improve'); 

 
 
 %% Output
 output = P; 
 
 






