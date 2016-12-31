function X=find_fprimary(A,b,iter_cg,epsilon)
[m,n]=size(b);
AT = zeros(m,n);
for i = 1:m
AT(i) = conj(A(i));	
end

S=[b];
X= rand(m,n);
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
function R=adjoint_operator(S,AT)
R = AT.*fft(S);



%%%%%%%%%#################%%%%%%%%%%%%
function Q=forward_operator(P,A)
Q = ifft(A.*P);