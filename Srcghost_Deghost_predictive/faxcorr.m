function output=faxcorr(x)

% x=x';

nx=length(x);
npad=floor(nx/2)+1;
padx=[zeros(1,npad) x zeros(1,npad)];


fpadx=fft(padx);
fpadoutput=conj(fpadx).*fpadx;


padoutput=real(ifft(fpadoutput));


output=padoutput(1:nx)';




