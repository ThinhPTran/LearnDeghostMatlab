function output = downward_extr(input, nf, omega, nkx, kx, nz, dz)

z = 0.0;
vwater = 1500.0;
F1 = zeros(size(input));
wavefield_fk = fft2(input);
nt = nf;
nx = nkx;

for i = 1:floor(nf/2)+1
  for j = 1:nkx
     if ((omega(i)*omega(i))/(vwater*vwater) - kx(j)*kx(j)>0)
        F1(i,j) = exp(1i*dz*sqrt((omega(i)*omega(i))/(vwater*vwater) - kx(j)*kx(j)));
     else
        F1(i,j) = 0.0; 
     end
  end
end


for i = 1:nz
   z = z+dz;
   wavefield_fk = wavefield_fk.*F1;
end


wavefield_fk(floor(nf/2)+1:nf,:) = 0.0;


for i = 2:floor(nf/2)+1
  wavefield_fk(nf+2-i,1) = conj(wavefield_fk(i,1));
end


for i = 2:floor(nf/2)+1
  for j = 2:nx
    wavefield_fk(nt+2-i,nx+2-j) = conj(wavefield_fk(i,j));
  end
end


output = real(ifft2(wavefield_fk));