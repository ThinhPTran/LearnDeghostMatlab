function out = qualmes(input)

finput = fft(input); 
atcor=real(ifft(finput.*conj(finput))); 
maxatcor = max(abs(atcor)); 

atcor = atcor/maxatcor; 

atcor(1)=0; 
%atcor(10:end)=0;

out = sum(abs(atcor)); 


