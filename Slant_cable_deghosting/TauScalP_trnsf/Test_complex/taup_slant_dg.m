function output = taup_slant_dg(recover_eps,alpha,vwater,nt,dt,npprime,dpprime,in) 


pdelta=2.0*sin(alpha)
npprimeshift=floor(pdelta/dpprime)


dtpprime=0.004;
nfpprime=npprime;
dfpprime=1.0/((nfpprime)*dtpprime);


fpprimemin=0.01;
fpprimemax=(nfpprime-2)*dfpprime;
res=4.0;
nfb=15;



t_lag=npprimeshift;
fil_lgth=30;


output=zeros(size(in));


for i_iter=1:nt


%% Test auto correlation
  tmpin=in(i_iter,:)';
  
  
  % Deconv
  tmpout=predctdeconv2(tmpin,fil_lgth,t_lag);
  
  
  % Rescale
%   maxin=max(abs(tmpin));
%   maxout=max(abs(tmpout));
%   ratio=maxin/maxout;
%   tmpout=tmpout*ratio;
  
  
  % Cleanup
  tmpout=cleanup(tmpin,tmpout,npprime,dtpprime,nfb,fpprimemin,fpprimemax,res);
  
  
  % Output ghost model
  tmpout=tmpin-tmpout;
  
  
  output(i_iter,:)=tmpout';
  
  
end








