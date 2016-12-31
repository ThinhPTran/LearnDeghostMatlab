function [output, freq_hz] = analyse_spec_fwd1(x,wmode,Nt,dt,nfb,fmin_hz,fmax_hz,Resolution)

Nh=floor(Nt/2);
df=1/(Nt*dt);
f=df*[0:Nh -Nh:-1]';


freq_hz=zeros(nfb,1);
fft_x=fft(x,Nt);


%% Prepare wavelet

% log interpolation
if(wmode==1) 
    fratio=fmin_hz/fmax_hz;
    for i_iter=1:nfb
        ww=(i_iter-1)/(nfb-1);
        freq_hz(i_iter)=fmax_hz*fratio^ww;
    end
end

% linear interpolation
if(wmode==2) 
    frange=fmax_hz-fmin_hz;
    for i_iter=1:nfb
        ww=(i_iter-1)/(nfb-1)*frange;
        freq_hz(i_iter)=fmax_hz-ww;
    end
end


sum_wavelet=zeros(Nt,1);
wavelet=zeros(Nt,nfb);


for i_iter=1:nfb
    for j_iter=1:Nh
        rarg=(freq_hz(i_iter)-f(j_iter))*Resolution/freq_hz(i_iter);
        if (rarg > 13)
            rarg=13;
        end
        wavelet(j_iter,i_iter)=exp(-0.5*rarg^2);
        sum_wavelet(j_iter)=sum_wavelet(j_iter)+wavelet(j_iter,i_iter);
    end
end


% Normalize the wavelet
for i_iter=1:nfb
    for j_iter=1:Nh
        wavelet(j_iter,i_iter)=wavelet(j_iter,i_iter)/sum_wavelet(j_iter);
    end
end



time_spectral=zeros(Nt,nfb);

for i_iter=1:nfb
    out_wavelet=zeros(Nt,1);
    for j_iter=1:Nh
        out_wavelet(j_iter,1)=wavelet(j_iter,i_iter)*fft_x(j_iter);
    end
    for j_iter=Nh+1:Nt
        out_wavelet(j_iter,1)=conj(out_wavelet(Nt-j_iter+2,1));
    end
    time_spectral(:,i_iter)=out_wavelet;
end

output = time_spectral;

