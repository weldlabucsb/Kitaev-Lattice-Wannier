function [J,U] =  numerics_testing(potentialDepth)
%% generate the lattice potential
%using the 6 cosine components of the lattice potential generated as output
%from Peter's code.
% clear all;
close all;
if (nargin < 1)

    potentialDepth = 10; %in Er
end

%     A = [1,1,0.6,0.5];
%     ph_deg = [0, 0, 90, -70];
%     th_deg = [0,90,180,270];
%     pol_deg = [0,0,0,0];
%     plots = 0; %boolean to turn plots on or off
    
%% Define the Lattice Paramters
%in terms of the lattice light wavelength (i.e. lambda = 1 here)
G = 4 * pi * [[1; 0] [0; 1]];
% The coordinates of the potential coefficients in the reciprocal lattice, given in units of
% the reciprocal lattice vectors, G
hkl = [[1; 0] [-1; 0] [0; 1] [0; -1]];
% The ratios of the corresponding potential coefficients
vi = [-1 -1 -1 -1];

%compute the real space fundamental lattice vectors:
A = compute_rsl_vectors(G); %in units of lattice light wavelength

%now using the n vector formalism to label the potential coordinates 
max_m = 5;
mLength = 2*max_m + 1;
Vcoeff = zeros(mLength,mLength);
for jj = 1:size(hkl,2)
    hklx = hkl(1,jj);
    hkly = hkl(2,jj);
    %conversion such that center is (0,0) in matlab indexing
    Vcoeff(hklx+(max_m+1),hkly+(max_m+1)) = Vcoeff(hklx+(max_m+1),hkly+(max_m+1)) + vi(jj);
end
rs_potential = fourier_to_real_fft(Vcoeff,max_m,50);
pot_max = max(rs_potential,[],'all');
pot_min = min(rs_potential,[],'all');
Vcoeff = Vcoeff*(potentialDepth/(pot_max-pot_min));


% toc
disp("%%%%%%%%%%%%%%% Constructing Hamiltonian %%%%%%%%%%%%%%%")
tic

%% Create Hamiltonian
%now this is where things start to get different from the other code, given
%that we are considering only certain values of the quasimomentum for
%computing the wannier states

L = 7; %how many real space lattice sites to use. MAKE AN ODD NUMBER
max_qm = floor(L/2); %to make sure that we are always in the first BZ
qsize = 2*max_qm + 1;
zone_number = 1; %how many zones to plot (for extended zone picture)
[quasiX,quasiY] = meshgrid((-max_qm:max_qm)./L,(-max_qm:max_qm)./L); %note that this should be considered more as the coefficient to the reciprocal lattice vector! not the actual magnitude itself
quasiX = quasiX.'; %due to the somewhat backwards meshgrid
quasiY = quasiY.'; %also these don't necessarily correspond to quasix and y anymore really. just reciprocal lattice vector coefficients. 
%hammy is the hamiltonian. Needs to accommodate all of the information in
%components but in 2-d
hammy = zeros(mLength^2,mLength^2,qsize,qsize);
for ii = 1:mLength % n1 (x)
    for jj = 1:mLength %n2 (y)
        for kk = 1:mLength % n1' (x')
            for ll = 1:mLength %n2' (y')
                for mm = 1:qsize
                    for nn = 1:qsize
                        if (ii==kk && jj==ll)
                            hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),mm,nn)=hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),mm,nn)+norm((quasiX(mm,nn)*G(:,1)+quasiY(mm,nn)*G(:,2)+(ii-(max_m+1))*G(:,1)+(jj-(max_m+1))*G(:,2))./(2*pi))^2;
                        end
                        if (abs(ii-kk) <= max_m && abs(jj-ll) <= max_m)
                              hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),mm,nn) = hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),mm,nn)+ Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
                        end
                    end
                end
            end
        end
    end
end

% for ii = 1:mLength
%     for jj = 1:mLength
%         for kk = 1:mLength
%             for ll = 1:mLength
%                 %to transform to the code, each of the indices actual
%                 %values in terms of the matlab indices are:
%                 % ii(actual) = ii - (max_m+1). This way the center of the
%                 % matrix corresponds to momentum index (0,0).
%                 if (ii==kk && jj==ll)
%                     hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:)=reshape(hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:),qsize,qsize)+((quasiX+ii-(max_m+1)).*(quasiX+ii-(max_m+1)))+((quasiY+jj-(max_m+1)).*(quasiY+jj-(max_m+1)));
%                 end
%                 if (abs(ii-kk) <= max_m && abs(jj-ll) <= max_m)
%                       hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:) = hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:)+ Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
%                 end
%             end
%         end
%     end
% end

num_bands = 4;
eigvals = zeros(mLength^2,mLength^2,qsize,qsize);
eigvecs = zeros(mLength^2,mLength^2,qsize,qsize);
toc
disp("%%%%%%%%%%%%Diagonalizing Hamiltonian %%%%%%%%%%%%%%%")
tic
for ii = 1:qsize
    for jj = 1:qsize
        %diagonalize for each point in the BZ
        [eigvecs(:,:,ii,jj),eigvals(:,:,ii,jj)] = eig(hammy(:,:,ii,jj));
    end
end
clear hammy;
eigvecs = eigvecs(:,1:num_bands,:,:);
eigvals = eigvals(1:num_bands,1:num_bands,:,:);
% keyboard;
%% De-compactify the eigenvectors:
%bands to be caluclated. Really there are only 2 non-trivial sublattice
%minima
bands = [1];
toc
disp("%%%%%%%%%%%%%%% De-Compactify Eigenvectors %%%%%%%%%%%%%%%%")
tic
%this is the eigenvector output from the eig function. Remember that we
%'compactified' the dimension here, so we want to make it back into the
%natural matrix form. This is really just for ease of use...
weightsMatrix = zeros(mLength,mLength,length(bands),qsize,qsize);
for kk = 1:length(bands)
    for ll = 1:qsize
        for mm = 1:qsize
            weightsMatrix(:,:,kk,ll,mm) = reshape(eigvecs(:,bands(kk),ll,mm),mLength,mLength).';
        end
    end
end
clear eigvecs;
%great! now we have a more interpretable to get the information that we
%need
%% Construct band projected position operators the analytic way...
toc 
disp("%%%%%%%%%%%% Construct Band-Projected Position Operators %%%%%%%%%%%%%%%%")
tic

m0 = ceil((L-1)./2)+0.5;

%note try to think about bsxfun here...
R1 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
R2 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
% ppm = ParforProgressbar(numel(R1), 'progressBarUpdatePeriod', 1.5);
parfor oo = 1:numel(R1)
    [row,col] = ind2sub([qsize.^2*length(bands),qsize.^2*length(bands)],oo);
    if (row >= col)
        [qy1,qx1,band1] = ind2sub([qsize,qsize,length(bands)],row);
        [qy2,qx2,band2] = ind2sub([qsize,qsize,length(bands)],col);
        [R1(oo),R2(oo)] = comp_elem(weightsMatrix,L,max_m,max_qm,band1,band2,qx1,qx2,qy1,qy2);
    end
%     ppm.increment();
end
% delete(ppm);
R1 = R1 + (R1-diag(diag(R1)))';
R2 = R2 + (R2-diag(diag(R2)))';

%% Diagonalize the Band Projected Position Operators.
toc
disp("%%%%%%%%%%%%%%%Diagonalize Band-Projected Position Operators%%%%%%%%%%%%%%%%")
tic
%First diagonalize 
[V2,D2] = eig(R2,'vector'); 
%now to sort
[D2,indices] = sort(D2,'comparisonMethod','real');
V2 = V2(:,indices);

%which x positions do we want to consider? Note these aren't the actual
%value of the real space x operator, but rather the inices of the latice
%unit cell. In this case, I want to take the one in the center as well as
%the one next to it so that I can compute tunelling element between them
%note that 1,1 is the bottom left position site
x_positions = [ceil(L/2) (ceil(L/2)+1)];
%and the same story for the y positions
y_positions = [ceil(L/2) (ceil(L/2)+1)];
%which wannier funcs do we want? For this lattice we have two, and you
%could choose either one. Make sure that they are valid indices
sub_unit_funcs = [1]; 
grouped_eigenvecs = zeros(length(D2),L,length(x_positions),length(sub_unit_funcs));
grouped_eigenvals = zeros(L,length(x_positions),length(sub_unit_funcs));


%now the new way to do things
for ii = 1:length(x_positions)
    %now within each position we need to group by sub-minima index wannier
    %function
    pos_eigvals = D2(1+L*length(bands)*(x_positions(ii)-1):L*length(bands)*x_positions(ii));
    pos_eigvecs = V2(:,1+L*length(bands)*(x_positions(ii)-1):L*length(bands)*x_positions(ii));
    for jj = 1:length(sub_unit_funcs)
        grouped_eigenvecs(:,:,ii,jj) = pos_eigvecs(:,1+L*(sub_unit_funcs(jj)-1):L*sub_unit_funcs(jj));
        grouped_eigenvals(:,ii,jj) = pos_eigvals(1+L*(sub_unit_funcs(jj)-1):L*sub_unit_funcs(jj));
    end
end

%partially localized states
% x_pos_index = 1;
% sub_index = 1;
% wannier_func_index = 1;
% wannier_func_realspace = bloch_to_real_fft(grouped_eigenvecs(:,sub_index,x_pos_index,wannier_func_index),U,qsize,bands);


%% Construct and Diagonalize Other Position Operator
toc
disp("%%%%%%%%%%%%%%%Diagonalize Subspace Position Operators%%%%%%%%%%%%%%%%")
tic
%Now we need to re-evaluate our basis and construct the other position
%operator in the new basis of the semi-degenerate subspace that we just
%created....

%note that the following can be vectorized! Just matrix multiply
%note that matlab does not allow for higher dimensional matrix
%multipliciation, although there is some code ndfun.c that we could try, I
%have that saved in case we want ultimate performance. This is already a
%vectorized calculation, and there isn't a straightforward way to deal with
%higher dimensional matrix multiplication in matlab
degen_operator = zeros(L,L,length(x_positions),length(sub_unit_funcs));
for ii = 1:length(x_positions)
    for jj = 1:length(sub_unit_funcs)
        degen_operator(:,:,ii,jj) = (grouped_eigenvecs(:,:,ii,jj)')*R1*(grouped_eigenvecs(:,:,ii,jj));
    end
end

wannier_states_in_bloch_basis = zeros(qsize.^2*length(bands),length(x_positions),length(y_positions),length(sub_unit_funcs));
wannier_states_positions = zeros(1,length(x_positions),length(y_positions),length(sub_unit_funcs));
%now go through and for each x position and sub_unit_func index, calculate
%each of the required y position states
for ii = 1:length(x_positions)
    for jj = 1:length(sub_unit_funcs)
        [V1,D1] = eig(degen_operator(:,:,ii,jj),'vector'); 
        [D1,indices] = sort(D1,'comparisonMethod','real'); %note there is no clustering here, and we don't have to bother with kmeans anymore
        V1 = V1(:,indices);
        for kk = 1:length(y_positions)
            %note here the states are in the basis of the grouped
            %eigenvectors that came out of the last caluclation. Therefore,
            %we have to put them back into the bloch basis!
            wannier_states_in_bloch_basis(:,ii,kk,jj) = grouped_eigenvecs(:,:,ii,jj)*V1(:,y_positions(kk));
            %we don't have to worry about any basis change for scalars!
            wannier_states_positions(1,ii,kk,jj) = D1(y_positions(kk));
        end
    end
end


%%Re-Phase the Bloch States
toc
disp("%%%%%%%%%%%%%%%%%%%%%%Re-Phase the Bloch states %%%%%%%%%%%%%%%%%%%")
tic
%here x_pos_ind and y_pos_ind are the wannier func in the center. 
x_pos_index = 1;
y_pos_index = 1;
wannier_func_index = 1;
wannier_func_realspace = bloch_to_real_fft(wannier_states_in_bloch_basis...
    (:,x_pos_index,y_pos_index,wannier_func_index),weightsMatrix,max_m,L,bands,0);
gamma = -angle(wannier_func_realspace(53,55));


for alpha = 1:length(bands)
    for qx = 1:qsize
        for qy = 1:qsize
            %note that in the language of what I have written down in the
            %notebook, this will be v^{\vec{c} = x_pos,y_pos}_{\vecq =
            %qx,qy, alpha = alpha}
            eigenvec_comp = wannier_states_in_bloch_basis((alpha-1)*(L)^2 + (qx-1)*L + qy,x_pos_index,y_pos_index,wannier_func_index);
            %now take k dot R
            k_dot_R = dot((x_positions(x_pos_index).*A(:,1) + y_positions(y_pos_index).*A(:,2)),...
                quasiX(qx,qy)*G(:,1)+quasiY(qx,qy)*G(:,2));
            weightsMatrix(:,:,alpha,qx,qy) = weightsMatrix(:,:,alpha,qx,qy).*exp(1i*(k_dot_R+gamma+angle(eigenvec_comp)));
        end
    end
end

% magnitudes_of_earlier_wannier =zeros(qsize.^2*length(bands),length(x_positions),length(y_positions),length(sub_unit_funcs));
% for xind = 1:length(x_positions)
%     for yind = 1:length(y_positions)
%         for subind = 1:length(sub_unit_funcs)
%             magnitudes_of_earlier_wannier(:,xind,yind,subind) = norm(wannier_states_in_bloch_basis(:,xind,yind,subind));
%         end
%     end
% end
%now that this is done, we can create invariable wannier functions anywhere
%we want on the lattice!
wannier_states_in_bloch_basis = construct_wannier_funcs(x_positions,y_positions,bands,sub_unit_funcs,qsize,A,G,quasiX,quasiY,L);

trying_change_mag = 0;
if (trying_change_mag)
    wannier_states_in_bloch_basis = wannier_states_in_bloch_basis.*magnitudes_of_earlier_wannier;
end
%%Generate Invariable Real-Space Plots
toc
disp("%%%%%%%%%%%%%%%%%%%%%%Generate Real Space Plots%%%%%%%%%%%%%%%%%%%")
tic

%now let's try to generate using the canonical definition with states with
%fixed phases.



if(0)
    wannier_plots(wannier_realspace_dft);
end

%% Calculate Bose Hubbard Parameters
toc
disp('%%%%%%%%%%%%%%%% Calculate Bose Hubbard Parameters %%%%%%%%%%%%%%%%%%%')
tic
%note that now we need to re-construct the hamiltonian in the sub-space of
%the bloch states that we are trying to consider. Since we already found
%the eigenbasis of the hamiltonian, this is just a question of construcing
%a diagonal matrix with the eigenvalues of the hamiltonian already found

%note that we only need the collection of bands that are specified in the
%calculation earlier.

sub_space_hamiltonian = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
vectorized_eigvals = zeros(qsize.^2,length(bands));
for ii = 1:length(bands)
   vectorized_eigvals(:,ii) = reshape(eigvals(bands(ii),bands(ii),:,:),qsize.^2,1);
end
sub_space_hamiltonian(:,:) = diag(reshape(vectorized_eigvals,(qsize.^2)*length(bands),1));




%calculate the J matrix element with the inner product of the hamiltonian
%and the wannier states in the bloch basis...
%in this case let's just calculate going one across in the x direction for
%the first wannier state
J = (wannier_states_in_bloch_basis(:,1,1,1)')*sub_space_hamiltonian*wannier_states_in_bloch_basis(:,2,1,1);

[wannier_realspace_dft,X,Y] = plotting_bloch(wannier_states_in_bloch_basis(:,1,1,1),weightsMatrix,max_m,L,bands,1,A);
x = X(1,:);
y = Y(:,1);
U = trapz(y,trapz(x,abs(wannier_realspace_dft).^4,2));
end