function y = convolution1(b,x)
[nb tmp] = size(b);
[nx tmp] = size(x);
ny = nb + nx - 1;
y = zeros(ny,1);
for ib=1:nb
    for ix=1:nx
        y(ib+ix-1) = y(ib+ix-1) + x(ix)*b(ib);
    end
end
y;