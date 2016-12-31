close all
clear all
clc

t = (0:0.004:8)';
f = (0:0.125:250)';
omega = 2*pi*f;

x = 1./cosh(100*(t-2));
plot(t,x);

b = zeros(2001,1);
b(1251) = 1;


figure(1);
plot(t,x);


%% Convolution (Method I)
conv1 = conv(b,x);

[m n] = size(conv1);
t11 = 0:0.004:(m-1)*0.004;

figure(2);
plot(t11,conv1);
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cross-correllation (Method I)
corr1 = xcorr(x,conv1);

[m n] = size(corr1);
t12 = 0:0.004:(m-1)*0.004;

figure(3);
plot(t12,corr1);
%%%%%%%%%%%%%%%%%%%%%%%%%


%% Convolution (Method II)
conv2 = convolution1(b,x);

[m n] = size(conv2);
t21 = 0:0.004:(m-1)*0.004;

figure(4);
plot(t21,conv2);
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cross-correllation (Method II)
corr2 = crosscorr1(x,conv2);

[m n] = size(corr2);
t22 = 0:0.004:(m-1)*0.004;

figure(5);
plot(t22,corr2);
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Convolution (Method III)
conv3 = convolution2(b,x);

[m n] = size(conv3);
t31 = 0:0.004:(m-1)*0.004;

figure(6);
plot(t31,conv3);
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cross-correllation (Method III)
corr3 = crosscorr2(b,conv3);

[m n] = size(corr3);
t32 = 0:0.004:(m-1)*0.004;

figure(7);
plot(t32,corr3);
%%%%%%%%%%%%%%%%%%%%%%%%%



