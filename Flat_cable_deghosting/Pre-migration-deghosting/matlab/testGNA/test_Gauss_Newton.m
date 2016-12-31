clear all
close all
clc

iter=1000;
epsilon = 0.01;
h = 0.01;


nt = 51
dt = 0.004
df = 1./((nt-1)*dt)
tmax = (nt-1)*dt
fmax = (nt-1)*df

t = (0:dt:tmax)';
f = (0:df:fmax)';
omega = 2*pi*f;


z = 20.0;
tau = 2.0*z/1500
index = floor(tau/dt) + floor(nt/2) + 1;

F = zeros(nt, 1);
F(floor(nt/2) + 1) = 1.0;
F(index) = -1.0;


%% Get primary
noise = randn(nt,1)/20.0;
wlet = ricker(40,0.004);
primary_sig =  zeros(nt,1);
for i = 1:13
   primary_sig(i+10) = wlet(i); 
end

figure(1)
plot(t,primary_sig);



%% Synthetize data with ghost
withghost = convolution2(primary_sig,F);% + noise;
figure(2)
plot(t,withghost);


%% Create Mirror data
mirror_data = -withghost;
figure(3)
plot(t,mirror_data);


fwithghost = fft(withghost);
figure(4);
plot(f,abs(fwithghost));


figure(5);
plot(f,angle(fwithghost));


correlation = xcorr(withghost,mirror_data);
[m n] = size(correlation);
tmaxtmp = (m-1)*dt;
t1 = (0:dt:tmaxtmp)';

figure(6);
plot(t1,correlation);



x = (1:nt)';
y= withghost;
betasize = nt + 2;

inita = 1.0;
inittau = 0.020;
initbeta = [primary_sig ; inita ; inittau];
J = zeros(nt,betasize);

bet = initbeta;
r = y - model(x,bet,dt);
err_new = r'*r;


for i_iter = 1:iter

    for i = 1:nt
        for j = 1:betasize
          tmpbet = bet;
          tmpbet(j) = bet(j) + h;
          J(i,j) = (model(x(i),tmpbet,dt) - model(x(i),bet,dt))/h;
        end
    end


    A = J'*J;
    b = J'*r;

    delta = conj_grad_solve(A,b,epsilon,1000);

    bet = bet + delta;
    
    
    r = y - model(x,bet,dt);
    err_old = err_new;
    err_new = r'*r;
    
    rms = sqrt(err_new);
    f = abs(err_new - err_old)/abs(err_old);
    
    if rms < epsilon ; break; end;
    if f < epsilon ; break; end;
end

i_iter
bet


y1 = model(x,bet,dt);









