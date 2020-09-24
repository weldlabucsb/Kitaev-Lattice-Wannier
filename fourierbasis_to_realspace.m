function [U] = fourierbasis_to_realspace(weights,max_m,xmat,ymat)
%NOTE: we have to be a bit careful when using meshgrid here. Unfortunately
%the A(i,j) matrix indexing convention conflicts with our traditional idea
%of (x,y) pairs. For (i,j), i represents a vertical displacement, but for (x,y)
%the first index represents horizontal displacement. If this isn't clear,
%try testing a simple example of meshgrid, or check MATLAB documentation
%site. Long story short we can just transpose the weights matrix to take
%this into account. That's why you see weights(jj,ii) not weights(ii,jj)
%also note that this function can be used for the bloch functions, since
%these are functions that are strictly periodic with the lattice...
[mxMat,myMat] = meshgrid(-max_m:max_m,-max_m:max_m);
U = zeros(length(xmat(:,1)),length(ymat(:,1)));
for ii = 1:(2*max_m+1)
    for jj = 1:(2*max_m+1)
        mx = mxMat(ii,jj);
        my = myMat(ii,jj);
        U = U + weights(jj,ii).*(exp(1i.*2.*pi.*mx.*xmat)).*(exp(1i.*2.*pi.*my.*ymat));
    end
end

end