%this is the meta code for running multiple iterations of the
%numerics_testing(and then eventually the lattice wannier code) to produce
%J plots vs the depth of the lattice 
n = 10;
depths = linspace(0.1,30,n);
Js = zeros(n,1);
Us = zeros(n,1);

for ii = 1:length(Js)
    [Js(ii),Us(ii)] = numerics_testing(depths(ii));
    disp([num2str(length(Js)-ii) ' runs remaining']);
end

hold on;
value = 2;
plot(depths,abs(Js),'linewidth',value);
set(gca,'YScale','log');


plot(depths,Us./Us(1),'linewidth',value);
set(gca,'YScale','log');

fontsize = 26;
xlabel('Lattice Depth [$E_{R}$]','interpreter','latex','fontsize',fontsize);

plot(depths,abs(Js)./(Us./Us(1)),'linewidth',value);
legend('J','U','J/U');


title(['Bose-Hubbard Model Parameters' ],'interpreter','latex','fontsize',fontsize);