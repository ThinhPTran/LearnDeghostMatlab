function y = convolution2(b,x)
[nb tmp] = size(b);
[nx tmp] = size(x);
ny = nx;
y = zeros(ny,1);
for iy=1:ny
    for ib=1:nb
        y(iy) = y(iy) + x(mod(iy - ib + floor(nb/2)+nx,nx)+1)*b(ib);
    end
end
y;