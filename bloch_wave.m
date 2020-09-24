function [U] = bloch_wave(eigvecs,max_m,xmat,ymat,quasiX,quasiY,bands)
%Now the goal of this function is to actually return the eigenfunctions of
%the Hamiltonian in real space. The other ones are the bloch functions (that have the
%periodicity of the lattice)
%This is the momenta grid.
[mxMat,myMat] = meshgrid(-max_m:max_m,-max_m:max_m);
mLength = 2*max_m + 1;
%Note that the first two indices are the real space representation of the
%bloch waves. The next two indices are the quasimomenta. The final one is
%the band index. I think that we should only need the first two, since we
%have two minima in a single lattice cell
U = zeros(length(xmat(:,1)),length(ymat(:,1)),length(bands),length(quasiX(:,1)),length(quasiY(:,1)));
for kk = bands
    %please see the note in the bloch function function about why this is
    %weightsMatrix(jj,ii) as opposed to (ii,jj).
    for ll = 1:length(quasiX(:,1))
        for mm = 1:length(quasiY(:,1))
            
            weightsColumn  = eigvecs(:,kk,ll,mm);
            weightsMatrix = zeros(mLength,mLength);
            %NOTE, this should be replaced with a reshape command!
            for ii = 1:mLength
                for jj = 1:mLength
                    weightsMatrix(ii,jj) = weightsColumn((mLength*(ii-1)+jj));
                end
            end
            
            for ii = 1:mLength
                for jj = 1:mLength
                    mx = mxMat(ii,jj);
                    my = myMat(ii,jj);
                    U(:,:,kk,ll,mm) = U(:,:,kk,ll,mm) +...
                        weightsMatrix(jj,ii).*(exp(1i.*2.*pi.*(mx+quasiX(ll,mm)).*xmat)).*(exp(1i.*2.*pi.*(my+quasiY(ll,mm)).*ymat));
                end
            end
        end
    end
    
end

end