clear all
close all
clc


%% Control parameters
vwater = 1500;
recover_eps = 0.05;
p_control = 1.0;
ctr_n = 1.0;
iter_cg = 1000;
res = 5.0;



%% Input data parameters
nx = 636;
nt = 2251;
ntau = nt;
np = 3001;

dt = 0.004;
dx = 12.5;
dtau = dt;
dp = 2.0/(vwater*(np+1));
fp = -1.0/vwater+dp;
fx = 100;


nf = nt;
nkx = nx;
df = 1.0/((nt)*dt);
dkx = 2.0*pi/((nx)*dx);


t=(0:dt:(nt-1)*dt)';
f=(0:df:(nf-1)*df)';
kx=dkx*[0:floor(nkx/2) -floor(nkx/2)+1:-1]';
x=(fx:dx:fx+(nx-1)*dx)';
tau=(0:dtau:(ntau-1)*dtau)';
p = (fp:dp:-fp)';
omega = 2.0*pi*f;


ntwindow=50;
wshift=40;
nw=fix((nt - ntwindow)/wshift) + 1;




%% Get orig and deghost
% fid=fopen('Thinh10760_in.bin','r');
% orig=fread(fid,[nt nx],'single');
% fclose(fid);
% 
% fid=fopen('Thinh10760_out.bin','r'); 
% deghost=fread(fid,[nt nx],'single'); 
% fclose(fid); 
% 
% tic
% tauporig=taup_fwd(nf,f,nx,x,np,p,orig); 
% taupdeghost=taup_fwd(nf,f,nx,x,np,p,deghost); 
% toc
% 
% 
% fid=fopen('Thinh10760_tauporig.bin','w'); 
% fwrite(fid,tauporig,'single'); 
% fclose(fid); 
% 
% fid=fopen('Thinh10760_taupdeghost.bin','w'); 
% fwrite(fid,taupdeghost,'single'); 
% fclose(fid); 


fid=fopen('Thinh10760_tauporig.bin','r'); 
tauporig=fread(fid,[nt np],'single'); 
fclose(fid); 

fid=fopen('Thinh10760_taupdeghost.bin','r'); 
taupdeghost=fread(fid,[nt np],'single'); 
fclose(fid); 


taupoutput=zeros(size(tauporig)); 

tic

for tridx = 1:np
  for twidx = 1:nw
        
    tridx
    twidx
        
    fidx = (twidx-1)*wshift + 1;
    lidx = (twidx-1)*wshift + ntwindow;

    tmporig = tauporig(fidx:lidx,tridx);
    tmpdeghost = taupdeghost(fidx:lidx,tridx); 
    
    
    mesin=qual_mes(tmporig); 
    mesgh=qual_mes(tmpdeghost); 
    
    
    if (mesgh>mesin)
      signout=sign(mean(tmpdeghost)); 
      
      tmpdeghost=(signout*tmporig+tmpdeghost)/2.0;  
    end
    
    
    for i_iter = 1:ntwindow
      weight = wshift./ntwindow;
      taupoutput(fidx + i_iter - 1,tridx)=taupoutput(fidx + i_iter - 1,tridx)+weight*tmpdeghost(i_iter);
    end
    

  end
end


toc


figure(); 
plot(1:nt,taupdeghost(:,50),1:nt,taupoutput(:,50)); 
legend('deghost','output'); 


fid=fopen('Thinh10760_taupclean.bin','w'); 
fwrite(fid,taupoutput,'single'); 
fclose(fid); 


