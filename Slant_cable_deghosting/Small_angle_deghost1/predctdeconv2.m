function output=predctdeconv2(in, fil_lgth, t_lag)

nt = length(in);
rin=faxcorr(in);

A=zeros(fil_lgth,fil_lgth);
for i_iter=1:fil_lgth
  A(i_iter:end,i_iter)=rin(1:fil_lgth-i_iter+1); 
end
A=A+tril(A,-1)';
b=rin(1+t_lag:fil_lgth+t_lag);


% Solve the problem
% disp('Solve the problem');
f=A\b;


f = [1,zeros(1,t_lag-1),-f']';
tmpout=conv(in,f);

output=tmpout(1:nt);








