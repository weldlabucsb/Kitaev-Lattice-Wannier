function [real_space,X,Y] = plotting_bloch(psi,weights_matrix,max_m,L,bands,sym,A)
%PLOTTING_BLOCH For plotting functions that are given in the bloch basis
%   This is very similar to the plotting_fourier function, except for bloch
%   functions. The point is to also provide X and Y matrices that allow for
%   accurate position labelling on the plots...

length1 = norm(A(:,1));
length2 = norm(A(:,2));
[X,Y] = meshgrid(linspace(-length1*L./2,length1*L./2,(2*max_m+1)*L),linspace(-length2*L./2,length2*L./2,(2*max_m+1)*L));

fourier_comps = zeros((2*max_m+1)*L);
%potentially bsxfun could be good here.
%https://www.mathworks.com/matlabcentral/answers/104127-how-to-multiply-multidimensional-arrays-with-a-column-vector
for kk = 1:length(bands)
    for ll = 1:L %quasimomentumx
        for mm = 1:L %quasimomentumy
            %here find the state in the fourier basis of the real space
            %lattice (since the fourier transform is linear, we can add the
            %vector before transforming it
            quasiMatrix = zeros(L);
            quasiMatrix(ll,mm) = 1;
            fourier_comps = fourier_comps + psi((kk-1)*(L)^2 + (ll-1)*L + mm).*kron(weights_matrix(:,:,kk,ll,mm),quasiMatrix);
        end
    end
end
if(sym)
    real_space = fftshift(ifft2(ifftshift(fourier_comps.'),'symmetric'));
else
    real_space = fftshift(ifft2(ifftshift(fourier_comps.')));
end
real_space = real_space.*(length(real_space).^2);
end

