function output = clean_non_physical_fk(ctr_n,vwater,nt,t,nf,f,omega,nx,x,np,p,nkx,kx,input) 

% Transform into FK domain
input_fk = fft2(input);


% Filtering
for i_iter = 1:floor(nf/2)+1
  for j_iter = 1:nkx
     if ((omega(i_iter)*omega(i_iter))/(ctr_n*vwater*vwater) - kx(j_iter)*kx(j_iter)<=0)
        input_fk(i_iter,j_iter) = 0;
     end
  end
end


for i_iter = 2:floor(nf/2)+1
  input_fk(nf+2-i_iter,1) = conj(input_fk(i_iter,1));
end


for i_iter = 2:floor(nf/2)+1
  for j_iter = 2:nkx
    input_fk(nf+2-i_iter,nkx+2-j_iter) = conj(input_fk(i_iter,j_iter));
  end
end


output = real(ifft2(input_fk));


