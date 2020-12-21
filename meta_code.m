%this is the meta code for running multiple iterations of the
%numerics_testing(and then eventually the lattice wannier code) to produce
%J plots vs the depth of the lattice 
n = 60;
depths = linspace(1,30,n);
Js = zeros(n,1);

for ii = 1:length(Js)
    Js(ii) = numerics_testing(depths(ii));
    disp([num2str(length(Js)-ii) ' runs remaining']);
end

plot(depths,abs(Js));
set(gca,'YScale','log');
xlabel('Lattice Depth, Er')
ylabel('One site tunneling element')
keyboard;