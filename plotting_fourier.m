function [U,X,Y] = plotting_fourier(weights,max_m,ext_m,A,times)
%PLOTTING_FOURIER Useful for plotting functions that are given in the
%fourier domain.
%   Does the same thing as fourier to real fft, however, also includes X
%   and Y matrices for a meshgrid plot such that things make sense with
%   units (in terms of the lattice light wavelength and all that)
%   NOTE right now this only really works for 2D. Will need to do some work
%   to extend this to higher dimensions if we want also note that I need to
%   think about extending this to the case where the reciprocal (and real)
%   lattice vectors aren't perpendicular. This might make things kind of
%   ugly, since you would have to take components and do some rouding.
%   Since this wouldn't necessarily make a lattice (i.e. if there are
%   potentially irrational angles) all the time. 
length1 = norm(A(:,1));
length2 = norm(A(:,2));
[X,Y] = meshgrid(linspace(-length1*times./2,length1*times./2,(2*ext_m+1)*times),linspace(-length2*times./2,length2*times./2,(2*ext_m+1)*times));
U = zeros(2*ext_m+1,2*ext_m+1);
U(ext_m-max_m+1:ext_m+max_m+1,ext_m-max_m+1:ext_m+max_m+1) = weights;
U = fftshift(ifft2(ifftshift(U.')));
U = U.*(length(U)^2);
U = repmat(U,times,times);
end

