function output = reshapespec(in, out, nx, nt, dt, nfb, fmin, fmax, res)

output = zeros(size(out));

sum_in=sum(in,2);
sum_out=sum(out,2);


%% Calculate gain parameters
[tspecin,tmp]=analyse_spec_fwd1(sum_in,1,nt,dt,nfb,fmin,fmax,res);
[tspecout,tmp]=analyse_spec_fwd1(sum_out,1,nt,dt,nfb,fmin,fmax,res);


ratio=zeros(1,nfb);

for i_iter=1:nfb
  maxin=max(abs(tspecin(:,i_iter)));
  maxout=max(abs(tspecout(:,i_iter)));
  ratio(i_iter) = maxout/maxin;
end


for i_iter = 1:nx
  
  tmp_out = out(:,i_iter);
  
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
  
  output(:,i_iter) = sum(tmp_tspecout,2);
  
end







