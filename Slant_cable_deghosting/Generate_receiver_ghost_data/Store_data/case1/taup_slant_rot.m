function output = taup_slant_rot(recover_eps,alpha,vwater,nt,dt,npprime,dpprime,in) 


pdelta=sin(alpha)


dtpprime=0.004;
nfpprime=npprime;
dfpprime=1.0/((nfpprime)*dtpprime);
fpprime=(0:dfpprime:(nfpprime-1)*dfpprime)';
omegapprime=2.0*pi*fpprime;


t_lag=pdelta*dtpprime/dpprime %+dtpprime


tshiftoper = exp(-1i*t_lag*omegapprime);


output=zeros(size(in));


for i_iter=1:nt


%% Test auto correlation
  tmpin=in(i_iter,:)';
  
  
  ftmpin=fft(tmpin);
  ftmpout=ftmpin.*tshiftoper;
  
  
  for j_iter=2:floor(nfpprime/2)
    ftmpout(nfpprime-j_iter+2)=conj(ftmpout(j_iter));    
  end
  
  tmpout=real(ifft(ftmpout));
  
  
  output(i_iter,:)=tmpout';
  
  
end








