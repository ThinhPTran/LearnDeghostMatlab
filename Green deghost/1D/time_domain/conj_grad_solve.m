function output=conj_grad_solve(ffilter,b,iter_cg)
[n, m]=size(b);
X= zeros(n,1);
R=[b];
Z=adjoint_operator(R,ffilter);
P=Z;
ZZ_new = Z'*Z;
err_old=R'*R;
for j = 1:iter_cg;
    Q=forward_operator(P,ffilter);
    alpha=ZZ_new/(Q'*Q);
    X=X+alpha*P;
    R=R-alpha*Q;
    Z=adjoint_operator(R,ffilter);
    ZZ_old=ZZ_new;
    ZZ_new=Z'*Z;
    beta=ZZ_new/ZZ_old;
    P=Z+beta*P;
    % error term computations
    err_new=R'*R;
    f = 200.0*(abs(err_old-err_new))/(err_old+err_new);
    rms=sqrt(err_new)
    err_old = err_new;
    if rms < 0.001 ; break; end;
    if f < 0.01; break;  end;
end

X=real(X); 
output=X;

num_iter=j


%%%%%%%%%#################%%%%%%%%%%%%
function Z=adjoint_operator(R,AT)
% adjffilter=transpose(ffilter)';
% Z = fft(adjffilter.*ifft(R));
% Z=transpose(ffilter)'.*R;
Z=AT*R;


%%%%%%%%%#################%%%%%%%%%%%%
function Q=forward_operator(P,A)
% Q = ifft(ffilter.*fft(P));
% Q=ffilter.*P;
Q=A*P;
