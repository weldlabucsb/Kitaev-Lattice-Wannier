function [wannier_blochspace_dft] = wannier_blochspace_construct(x_pos_index,y_pos_index,A,G,bands,qsize,x_positions,y_positions,quasiX,quasiY,L)
%WANNIER_BLOCHSPACE_CONSTRUCT Construct the wannier states in the bloch
%basis for the specified location
%   Detailed explanation goes here
wannier_blochspace_dft = zeros(qsize.^2*length(bands),1);

for alpha = 1:length(bands)
    for qx = 1:qsize
        for qy = 1:qsize
            k_dot_R = dot((x_positions(x_pos_index).*A(:,1) + y_positions(y_pos_index).*A(:,2)),...
                quasiX(qx,qy)*G(:,1)+quasiY(qx,qy)*G(:,2));
            wannier_blochspace_dft((alpha-1)*(L)^2 + (qx-1)*L + qy) = exp(-1i*k_dot_R);
        end
    end
end


wannier_blochspace_dft = wannier_blochspace_dft./norm(wannier_blochspace_dft);

end

