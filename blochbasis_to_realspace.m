function [real_space] = blochbasis_to_realspace(vector,U,qsize,bands)
%for now this is just for one state, and will take you from the bloch basis
%to real space
%the input U is a collection of the real space bloch waves, generated
%using the bloch_wave function below
%find the dimensions of the real space representation of the collection of
%bloch states (the first two dimensions of the matrix U)
dims = size(U);
real_space = zeros(dims(1:2));
%potentially bsxfun could be good here.
%https://www.mathworks.com/matlabcentral/answers/104127-how-to-multiply-multidimensional-arrays-with-a-column-vector
for ii = bands
    for jj = 1:qsize
        for kk = 1:qsize
            real_space(:,:) = real_space(:,:) +  vector((ii-1)*(qsize)^2 + (jj-1)*qsize + kk)...
                .*U(:,:,ii,jj,kk);
        end
    end
end

end