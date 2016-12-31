function output = slant_tauscp_bwd(iter_cg,vwater,x0,x_noff,alpha,nf,f,nx,x,np,p,input) 


finput = fft(input);

A = zeros(np,nx);
output = zeros(nf,np);
foutput = zeros(nf,nx);
tmp_in = zeros(np,1);
tmp_out = zeros(nx,1);



%% Inverse for original data
for i_iter = 1:floor(nf/2)
  
  tmp_in = transpose(finput(i_iter,:));
    
  for j_iter = 1:np
    for k_iter = 1:nx
      A(j_iter,k_iter) =...
      exp(1i*2.0*pi*f(i_iter) ...
          *(p(j_iter)*x(k_iter) ... 
            + 2*(x(k_iter)+x0) ...
               *tan(alpha) ...
               *sqrt(1-p(j_iter)*p(j_iter)*vwater*vwater) ...
               /vwater ...
            - 2*(x_noff+x0)*tan(alpha)/vwater ...
           ) ...
         ); 
    end
  end
  
  tmp_out = conj_grad_solve(A,tmp_in,iter_cg);
   
  foutput(i_iter,:) = transpose(tmp_out);
  
end


for i_iter = 2:floor(nf/2)
  foutput(nf-i_iter+2,:) = conj(foutput(i_iter,:)); 
end


output = real(ifft(foutput));



