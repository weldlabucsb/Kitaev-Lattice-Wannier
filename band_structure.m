%now to do the lattice bands!!


lattice_wavelength = 1.064e-6;%meters
num_band = 5;
lattice_depth = 20; %Er
v_0 = lattice_depth;
max_m = 50;
m_len = max_m*2+1

%here q is the quasimomentum, but perhaps I should get into the habit of
%making this k since I believe that this is the more standard notations
qs = linspace(-1,1,100);
%initialize the matrix
h = zeros(max_m*2+1,max_m*2+1,length(qs));
for ii = 1:m_len
    for jj = 1:m_len
        if ii == jj
            %these are the diagonal terms
            h(ii,jj,:) = (v_0)./2 + (qs+2*(max_m+1-ii)).^2;%here it is actually q/(light wavevector), so brilluoin zone is -1 to 1
        end
        if ii == jj+1 || ii == jj-1
            %the off-diagonal 'coupling' terms. This makes the system
            %non-trivial
            h(ii,jj,:) = (v_0)./4;
        end
    end
end
%initialize the array to store the eigenvalues
eigs = zeros(m_len,length(qs));
for i = 1:length(qs)
    %find the actual band sturucture by compting the energies.
    eigs(:,i) = eig(h(:,:,i));
end
figure;
hold all;
for i = 1:num_band
    %plot the energy vs quasimomentum
    plot(qs,eigs(i,:));
end
xlabel('quasimomentum(q), [k_{l}]')
ylabel('energy in Recoils')

%now we want to be able to plot the bloch functions... We need the
%coefficients of the eigenvector to do this since this is how we get the
%function...

%since we need to do this at a specific quasi-momentum, let us do it at
%zero for now
q_index = 50;
q = qs(q_index)
%to find the bloch functions we actually need the eigenvectors, not just
%the eigenvalues
[weights,eigs] = eig(h(:,:,q_index));
m_vec = -max_m:max_m;
x = linspace(0,20,200);
u = zeros(1,length(x));
for ii = 1:length(x)
    u(ii) = abs(bloch_func(weights(:,1),m_vec,x(ii)));
end
figure;
plot(x,u);

function [u] = bloch_func(weights,m_vec,x)%x is in units of d, the lattice spacing (pi/kl)
u = 0;
    for ii = 1:length(m_vec)
        m = m_vec(ii);
        u = u+weights(ii).*(exp(1i.*.2.*m.*pi.*x));
    end
end



%now to try my way. The above code works

% 
% function eigs = find_eigs(q,max_m,v_0)
% ms = linspace(-max_m,max_m);
% h = zeros(length(ms));
% for ii = 1:length(ms)
%     for jj = 1:length(ms)
%         if ii == jj
%             h(ii,jj) = (v_0)./2 + (q+2*(ms(ii))).^2;%here it is actually q/(lattice constant), so brilluoin zone is -1 to 1
%         end
%         if ii == jj+1 || ii == jj-1
%             h(ii,jj) = (v_0)./4;
%         end
%     end
% end
% eigs = eig(h);
% end

