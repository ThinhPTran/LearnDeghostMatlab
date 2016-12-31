function y = model(x,beta,dt)

[n m] = size(x);
y = zeros(n,1);

[tmpnt m] = size(beta);

nt = tmpnt -2;

a = beta(nt+1)
taushift = fix(beta(nt+2)/dt)

primary = zeros(nt,1);
for i = 1:nt
    primary(i) = beta(i);
end

ghost = -circshift(primary,[taushift 0]);
for i = 1:taushift
   ghost(i) = 0; 
end

withghost = primary+a*ghost;

for i = 1:n
   y(i) = withghost(x(i)); 
end

