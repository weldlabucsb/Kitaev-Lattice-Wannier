%find the bloch_functions



% figure;
% max_m = 20;
% lattice_depth = 0;
% n = 5;
% q = linspace(-1,1,100);
% eigs = zeros(n,length(q));
% for ii = 1:length(q)
%     [vectors,vals] = optical_potential(q(ii),max_m,lattice_depth);
%     eigs(:,ii) = vals(1:n);
% end
% hold all;
% for ii = 1:n
%     plot(q,eigs(ii,:));
% end

%that was just calculating the band strcuture
figure(23);
max_m = 600;
m_vec = -max_m:max_m;
lattice_depth = 100;
n = 1; %the lowest band (because matlab)
x = linspace(0,1,1000);
u = zeros(1,length(x));
[vectors,vals] = optical_potential(0,max_m,lattice_depth);
for ii = 1:length(x)
    u(ii) = real(bloch_func(vectors(:,n),m_vec,x(ii)));
end
plot(x,u);
xlabel('lattice sites')
ylabel('Bloch Function u')
title('Bloch fxn in position space')



function [u] = bloch_func(weights,m_vec,x)%x is in units of d, the lattice spacing (pi/kl)
    u = 0;
    for ii = 1:length(m_vec)
        m = m_vec(ii);
        weights(ii);
        arg = 1j*2*m*pi*x;
        exponed = exp(arg);
        thing = weights(ii).*exponed;
%         u = thing
        u = u+weights(ii).*exponed;
%         u=u*conj(u);
%         keyboard
%         u = 
%why do the above lines (specifically exp(arg)) always return an exponed
%thing of one??
%         m = m_vec(ii);
%         disp('x')
%         x
%         arg = 2*m*pi*x
%         disp('exponed')
%         exponed = cos(arg) + 1j*sin(arg)
%         u = u + weights(ii).*exponed;
    end
    u=conj(u).*u;
end

%ok, so now I have a function that gives me the eigenvalues as a column
%vector and the components of the complex fourier basis as a matrix
function [vectors,vals] = optical_potential(q,max_m,lattice_depth)
lattice_wavelength = 1.064e-6;%meters
num_band = 5;
v_0 = lattice_depth;
m_len = max_m*2+1
m_vec = -max_m:max_m;

h = zeros(m_len,m_len);
for ii = 1:m_len
    for jj = 1:m_len
        if ii == jj
            h(ii,jj) = (v_0)./2 + (q+2*(m_vec(ii))).^2;%here it is actually q/(lattice constant), so brillouin zone is -1 to 1
        end
        if ii == jj+1 || ii == jj-1
            h(ii,jj) = (v_0)./4;
        end
    end
end
[vectors,vals] = eig(h);%vectors is a matrix while vals is a column vector
vals = eig(h);
end