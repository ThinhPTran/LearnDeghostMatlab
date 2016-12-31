close all

x = (0:0.001:0.02);
% hist(storestatistic);
[n xout] =  hist(storestatistic,x);
normn = n./sum(n);
bar(xout,normn,'histc');
xlabel('$tau$','interpreter','latex','fontsize',18);
ylabel('$Percentage$','interpreter','latex','fontsize',18);
title('Normalized histogram for found tau in two objectives case Namp/Sigamp = 0.05');
set(gca,'fontsize',14); % increase font size
saveas(gcf,'histogram_for_found_tau_in_two_objectives_case_0_05.eps','eps2c');