function output = qual_mes(input)

finput=fft(input); 

autocor=real(ifft(transpose(finput').*finput)); 
% dautocor=diff(autocor); 
dautocor=diff(diff(autocor)); 
maxdau=max(abs(dautocor)); 
dautocor=dautocor/maxdau; 


nt=length(dautocor);
nh=floor(nt/2); 

% augdautocor=dautocor(1:nh).*(1:nh).*(1:nh);
augdautocor=dautocor(1:nh); 

output=sum(abs(augdautocor)); 











