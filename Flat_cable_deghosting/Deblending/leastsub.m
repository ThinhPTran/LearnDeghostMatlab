function output = leastsub(input1,input2, nx, nt, epsilon, iter_cg)

tmp_b = zeros(nt,1);
tmp_A = zeros(nt,nt);
output = zeros(nt,1);

for i = 1:nx
   i
   tmp_b(:,1) = input1(:,i);
   tmp_A = diag(input2(:,i));
   A = tmp_A + epsilon*eye(nt);
   b = tmp_b;
   x = conj_grad_solve(A,b,iter_cg);
   tmp_output = b-A*x;
   output(:,i) = tmp_output(:);
end
