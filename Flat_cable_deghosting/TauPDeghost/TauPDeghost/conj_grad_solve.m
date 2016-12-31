function X=conj_grad_solve(A,b,iter_cg)
[m,n]=size(A);
AT = A';
X= zeros(n,1);
S=[b];
R=adjoint_operator(S,AT);
P=R;
rr_new = R'*R;
err_old=S'*S;
for j = 1:iter_cg;
    Q=forward_operator(P,A);
    alpha=rr_new/(Q'*Q);
    X=X+alpha*P;
    S=S-alpha*Q;
    R=adjoint_operator(S,AT);
    rr_old=rr_new;
    rr_new=R'*R;
    beta=rr_new/rr_old;
    P=R+beta*P;
    % error term computations
    err_new=S'*S;
    f = 200.0*(abs(err_old-err_new))/(err_old+err_new);
    rms=sqrt(err_new);
    err_old = err_new;
    if rms < 0.001; break; end;
    if f < 0.01; break;  end;
end

% j
% rms
% f

%%%%%%%%%#################%%%%%%%%%%%%
function R=adjoint_operator(S,AT)
R = AT*S;


%%%%%%%%%#################%%%%%%%%%%%%
function Q=forward_operator(P,A)
Q = A*P;
