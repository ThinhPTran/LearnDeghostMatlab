close all
clear all
clc

t = [0:0.004:8]';
f = [0:0.125:250]';



fid = fopen('primary.bin','r');
primary = fread(fid,[2001 1001],'single');
fclose(fid);

fid = fopen('prighost.bin','r');
prighost = fread(fid,[2001 1001],'single');
fclose(fid);

fid = fopen('prighostout.bin','r');
prighostout = fread(fid,[2001,1001],'single');
fclose(fid);


figure(1)
plot(t,prighost(:,500),'blue',t,prighostout(:,500),'red');
legend('input data','deghosted data');


h = spectrum.welch('hamming',64);
h.WindowName='Chebyshev';
hpsd1 = psd(h,prighost(:,500),'Fs',250);
Pxx1 = hpsd1.Data;
W = hpsd1.Frequencies;

hpsd2 = psd(h,prighostout(:,500),'Fs',250);
Pxx2 = hpsd2.Data;

hpsd = dspdata.psd([Pxx1, Pxx2],W,'Fs',250)


figure(2)
plot(hpsd);
legend('input data','deghosted data')


absfinput = zeros(2001,1);
absfdeghosted = zeros(2001,1);

for i = 1:1001
absfinput = absfinput +  abs(fft(prighost(:,i)));
absfdeghosted = absfdeghosted + abs(fft(prighostout(:,i)));
end

absfinput = absfinput./1001;
absfdeghosted = absfdeghosted./1001;

figure(3)
plot(f,absfinput,'blue',f,absfdeghosted,'red');
legend('spec input','spec deghosted');



