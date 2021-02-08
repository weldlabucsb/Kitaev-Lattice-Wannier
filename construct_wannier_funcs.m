function [wannier_states_in_bloch_basis] = construct_wannier_funcs(x_positions,y_positions,bands,sub_unit_funcs,qsize,A,G,quasiX,quasiY,L)
%CONSTRUCT_WANNIER_FUNCS For constructing the wannier functions once the
%Bloch states are re-phased
%   Ouputs a matrix of vectors which are wannier states in the bloch basis.
%   

wannier_states_in_bloch_basis = zeros(qsize.^2*length(bands),length(x_positions),length(y_positions),length(sub_unit_funcs));

for ii = 1:length(x_positions)
    for jj = 1:length(y_positions)
        for cell_ind = 1:length(sub_unit_funcs)
            for alpha = 1:length(bands)
                for qx = 1:qsize
                    for qy = 1:qsize
                        k_dot_R = dot((x_positions(ii).*A(:,1) + y_positions(jj).*A(:,2)),...
                            quasiX(qx,qy)*G(:,1)+quasiY(qx,qy)*G(:,2));
                        %note I realize that the right phases for the bloch
                        %states DO depend on the function of our choosing (i.e. the band and the sub-unit cell index).
                        %This needs to be done relatively carefully, but
                        %for now is fine I think
                        wannier_states_in_bloch_basis((alpha-1)*(L)^2 + (qx-1)*L + qy,ii,jj,cell_ind)=exp(-1i*k_dot_R);
                    end
                end
            end
        end
    end
end


end

