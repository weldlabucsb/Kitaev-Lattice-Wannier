function [real_space] = bloch_to_real_fft(psi,weights_matrix,max_m,L,bands,sym)
%using the tensor product method, go from the bloch basis to realspace.
%This will require the set of bloch functions and the associated
%quasimomentum mesh. 
fourier_comps = zeros((2*max_m+1)*L);
%potentially bsxfun could be good here.
%https://www.mathworks.com/matlabcentral/answers/104127-how-to-multiply-multidimensional-arrays-with-a-column-vector
for kk = 1:length(bands)
    for ll = 1:L
        for mm = 1:L
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