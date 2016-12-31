function output = clean_non_physical_fk(ctr_n,vwater,nt,t,nf,f,omega,nx,x,np,p,nkx,kx,input1,input2) 

% Transform into FK domain
input1_fk = fft2(input1);
input2_fk = fft2(input2);


% Filtering
for i_iter = 1:floor(nf/2)+1
  for j_iter = 1:nkx
     if ((omega(i_iter)*omega(i_iter))/(ctr_n*vwater*vwater) - kx(j_iter)*kx(j_iter)<=0)
        input1_fk(i_iter,j_iter) = 0;%input2_fk(i_iter,j_iter);
     end
  end
end


for i_iter = 2:floor(nf/2)+1
  input1_fk(nf+2-i_iter,1) = conj(input1_fk(i_iter,1));
end


for i_iter = 2:floor(nf/2)+1
  for j_iter = 2:nkx
    input1_fk(nf+2-i_iter,nkx+2-j_iter) = conj(input1_fk(i_iter,j_iter));
  end
end


output = real(ifft2(input1_fk));


