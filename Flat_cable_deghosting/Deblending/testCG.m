

iter_cg = 100;

A=[ 2 -1 4 6i; 3 1 1 7; 9 -1 -2 -1]
b = [ 1 ; 4i ; -3]


X=conj_grad_solve(A,b,iter_cg)

