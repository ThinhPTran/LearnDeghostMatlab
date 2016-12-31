function y = gparabol(x,beta)

[m n] = size(x);
myone = ones(m,1);

y = beta(1)*x.*x + beta(2)*x + beta(3)*myone;