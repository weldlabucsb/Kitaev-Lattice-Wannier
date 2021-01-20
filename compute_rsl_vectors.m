function [R] = compute_rsl_vectors(G)
%compute_rsl_vectors Compute the real space lattice basis vectors from the
%reciprocal lattice basis vectors
if (length(G) == 2)
    %if in 2-D, have to add fake 3D to make cross products feasible
    G = [G;[0,0]];
    G = [G,[0;0;1]];
    
    R = zeros(2);
    c1 = 2.*pi.*cross(G(:,2),G(:,3))./det(G);
    R(:,1) = c1(1:2);
    c2 = 2.*pi.*cross(G(:,3),G(:,1))./det(G);
    R(:,2) = c2(1:2);
elseif (length(G) ==3)
    R = zeros(3);
    R(:,1) = 2.*pi.*cross(G(:,2),G(:,3))./det(G);
    R(:,2) = 2.*pi.*cross(G(:,3),G(:,1))./det(G);
    R(:,3) = 2.*pi.*cross(G(:,1),G(:,2))./det(G);
end

end

