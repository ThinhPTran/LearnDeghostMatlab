close all
clear all
clc


eps = 0.1;
vwater = 1500;
amplitude = 1;


nt = 501;
n2 = 240;
dt = 0.004;

df = 1./((nt-1)*dt);
tmax = (nt-1)*dt;
fmax = (nt-1)*df;

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;



%% t for checking
z = zeros(n2,1); 
zout = zeros(n2,1); 
tau = zeros(n2,1);


for iter = 1:n2
  z(iter) = (iter-1)*0.3125+5.0;  
  tau(iter) = 2.0*z(iter)/vwater; 
end


%% Get data
noise = 0.1*amplitude*randn(nt,1);


fid = fopen('slantcable_03125.bin','rb');
datain = fread(fid,[nt n2],'single');
fclose(fid);

for iter = 1:n2
   datain(:,iter) = datain(:,iter) + noise;  
end


P = zeros(size(datain));


testiter=100

%% Deghost 
for iter = 1:n2
iter
[P(:,iter) zout(iter)] = deghostfunc(datain(:,iter), z(iter), vwater, eps, nt, dt); 
end


zout;

 
figure();
plot(t,datain(:,testiter),'blue',t,P(:,testiter),'red');
legend('withghost','Primary recovered');


figure();
plot(1:n2,z,1:n2,zout);
legend('z','zout'); 


finput=abs(fft(datain(:,testiter)));
foutput=abs(fft(P(:,testiter)));

finputdb = mag2db(finput);
foutputdb = mag2db(foutput);


figure();
plot(f,finput,f,foutput);
legend('finput','foutput'); 


figure();
plot(f,finputdb,f,foutputdb); 
legend('finputdb','foutputdb'); 


fid = fopen('slantcable_03125_in.bin','wb');
fwrite(fid,datain,'single');
fclose(fid);


fid = fopen('slantcable_03125_out.bin','wb');
fwrite(fid,P,'single');
fclose(fid);















