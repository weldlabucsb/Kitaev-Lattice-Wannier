function [U] = fourier_to_real_fft(weights,max_m,ext_m)
%Very similar to fourierbasis_to_realspace, but using FFT to make things
%faster and more accurate. However, the tradeoff is that we are only able
%to see one unit cell here. 
%max_m is the highest momentum state that is in weights, and the fourier
%components before the transform are padded up to ext_m with zeros.
%This is to make the output plot smoother without changing the actual plot.
U = zeros(2*ext_m+1,2*ext_m+1);
U(ext_m-max_m+1:ext_m+max_m+1,ext_m-max_m+1:ext_m+max_m+1) = weights;
U = ifft2(ifftshift(U.'));
U = U.*(length(U)^2);
end