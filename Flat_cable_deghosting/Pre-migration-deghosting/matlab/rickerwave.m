function P=rickerwave(fc,tc,t);
P=(1-2*pi^2*(fc)^2*(t-tc).^2).*exp(-pi^2*(fc)^2.*(t-tc).^2);
end