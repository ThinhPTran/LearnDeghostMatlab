function X = forward_PEF(A,iter_cg,epsilon)
% Filter of form [f1 f2]
[m,n]=size(A);
S=circshift(A,[-1 0]);
X= zeros(2,1);
R=adjoint_operator(S,A,m);
P=R;
rr_new = R'*R;
err_old=S'*S;
for j = 1:iter_cg;
    Q=forward_operator(P,A,m);
    alpha=rr_new/(Q'*Q);
    X=X+alpha*P;
    S=S-alpha*Q;
    R=adjoint_operator(S,A,m);
    rr_old=rr_new;
    rr_new=R'*R;
    beta=rr_new/rr_old;
    P=R+beta*P;
    %% error term computations
    err_new=S'*S;
    f = 200.0*(abs(err_old-err_new))/(err_old+err_new)
    rms=sqrt(err_new)
    err_old = err_new;
    if rms < epsilon ; break; end;
    if f < epsilon; break;  end;
end
j

%%%%%%%%%#################%%%%%%%%%%%%
function R=adjoint_operator(S,A,m)
tmp = crosscorr2(A,S);
R = zeros(2,1);
R(1) = tmp(floor(m/2) + 1);
R(2) = tmp(floor(m/2) + 2);




%%%%%%%%%#################%%%%%%%%%%%%
function Q=forward_operator(P,A,m)
tmp = zeros(m,1);
tmp(floor(m/2) + 1) = P(1);
tmp(floor(m/2) + 2) = P(2);
Q = convolution2(tmp,A);
