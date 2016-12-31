function output = taup_slant_dg(recover_eps,alpha,vwater,fil_lgth,nt,dt,np,dp,in) 


fp=-1.0/vwater+dp;
p=(fp:dp:-fp)';
pdf=1.0/(np*dp);
pf=(0:pdf:(np-1)*pdf);
pomega=2.0*pi*pf;


pdelta=2.0*sin(alpha)/(vwater) 
npshift=pdelta/dp

t_lag=npshift;
% t_lag=87;


output=zeros(size(in));


axcorr=zeros(nt,1);


for i_iter=1:nt
% i_iter=125;


%% Test auto correlation
  data=in(i_iter,:)';
  tmp_out=faxcorr(data);
  
  
  [val ind]=max(abs(tmp_out(50:250)));
  axcorr(i_iter)=ind+49;
  
  
  t_lag=axcorr(i_iter);
  

%% Take one time slice
  tmpout = predctdeconv2(data,fil_lgth,t_lag);

  
  output(i_iter,:)=tmpout';
  
  
%   t_lag
%   figure();
%   plot(p,data,p,tmpout);
%   legend('data','tmpout');
  
  
end


figure();
plot(axcorr);
legend('axcorr');






