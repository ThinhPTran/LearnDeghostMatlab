function output = taup_slant_dg(recover_eps,alpha,vwater,nt,dt,np,dp,in) 

pdf=1.0/(np*dp);
pf=(0:pdf:(np-1)*pdf);
pomega=2.0*pi*pf;


pdelta=sin(2.0*alpha)/(vwater) 
npshift=pdelta/dp
% pdelta=100*dp;



output=zeros(size(in));

for i_iter=1:nt
   
    data=in(i_iter,:);
    fdata=fft(data);
    fout=zeros(size(fdata));
    
    minustshift=exp(1i*pdelta*pomega);
    F2=1-minustshift;
    denominator=2.0-2.0*cos(pomega*pdelta)+recover_eps;
    
    fout =(fdata.*F2)./denominator;
    
    for j_iter = np:-1:floor(np/2)
        fout(j_iter) = conj(fout(np-j_iter+2));
    end
    
    output(i_iter,:)=real(ifft(fout));
    
end







