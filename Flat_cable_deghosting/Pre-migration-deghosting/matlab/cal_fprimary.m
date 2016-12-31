function X=cal_fprimary(A,W,b,iter_cg,epsilon)
[m,n]=size(b);
AT = zeros(m,n);
WT = W;
for i = 1:2001
AT(i) = conj(A(i));	
endfor
S=[b];
X= zeros(m,n);
R=adjoint_operator(S,AT,WT);
P=R;
rr_new = R'*R;
err_old=S'*S;
for j = 1:iter_cg;
    Q=forward_operator(P,A,W);
    alpha=rr_new/(Q'*Q);
    X=X+alpha*P;
    S=S-alpha*Q;
    R=adjoint_operator(S,AT,WT);
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
function R=adjoint_operator(S,AT,WT)
R = WT.*AT.*S;



%%%%%%%%%#################%%%%%%%%%%%%
function Q=forward_operator(P,A,W)
Q = W.*A.*P;
