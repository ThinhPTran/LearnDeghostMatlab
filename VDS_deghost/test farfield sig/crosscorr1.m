function b = crosscorr1(x,y)
[ny tmp] = size(y);
[nx tmp] = size(x);
nb = ny - nx + 1;
b = zeros(nb,1);
for ib=1:nb
    for ix=1:nx
        b(ib) = b(ib) + x(ix)*y(ib+ix-1);
    end
end
b;