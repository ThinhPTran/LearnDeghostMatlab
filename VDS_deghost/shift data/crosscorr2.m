function b = crosscorr2(x,y)
[ny tmp] = size(y);
[nx tmp] = size(x);
nb = ny;
b = zeros(nb,1);
for iy=1:ny
    for ib=1:nb
        b(ib) = b(ib) + x(mod(iy - ib + floor(nb/2)+nx,nx)+1)*y(iy);
    end
end
b;