function output = reshapespec(in, out, nt, dt, nfb, fmin, fmax, res)

output = zeros(size(out));


%% Calculate gain parameters
[tspecin,tmp]=analyse_spec_fwd1(in,1,nt,dt,nfb,fmin,fmax,res);
[tspecout,tmp]=analyse_spec_fwd1(out,1,nt,dt,nfb,fmin,fmax,res);


ratio=zeros(1,nfb);

for i_iter=1:nfb
  maxin=max(abs(tspecin(:,i_iter)));
  maxout=max(abs(tspecout(:,i_iter)));
  ratio(i_iter) = maxout/maxin;
end


  
tmp_out = out;
  
[tspecout,tmp] = analyse_spec_fwd1(tmp_out,1,nt,dt,nfb,fmin,fmax,res);
    
for j_iter = 1:nfb
  if (ratio(j_iter)>0)
    tspecout(:,j_iter)=tspecout(:,j_iter)/ratio(j_iter); 
  end
end
  
tmp_tspecout=zeros(size(tspecout));

for j_iter = 1:nfb
  tmp_tspecout(:,j_iter)=real(ifft(tspecout(:,j_iter)));
end
  
output=sum(tmp_tspecout,2);
  








