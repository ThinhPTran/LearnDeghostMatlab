function output = taup_fwd(nf,f,nx,x,np,p,input) 

finput = fft(input);

A = zeros(np,nx);
output = zeros(nf,np);
foutput = zeros(nf,np);
tmp_in = zeros(nx,1);
tmp_out = zeros(np,1);


for i_iter = 1:floor(nf/2)
  
  tmp_in = transpose(finput(i_iter,:));
    
  for j_iter = 1:np
    for k_iter = 1:nx
      A(j_iter,k_iter) = exp(1i*2.0*pi*f(i_iter)*x(k_iter)*p(j_iter));   
    end
  end
  
  tmp_out = A*tmp_in;
  
  foutput(i_iter,:) = transpose(tmp_out);
  
end


for i_iter = 2:floor(nf/2)
  foutput(nf-i_iter+2,:) = conj(foutput(i_iter,:)); 
end


output = real(ifft(foutput));










