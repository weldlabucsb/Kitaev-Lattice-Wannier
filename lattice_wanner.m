function[] = lattice_wannier()
%% generate the lattice potential
%using the 6 cosine components of the lattice potential generated as output
%from Peter's code.
A = [1,1,0.6,0.5];
ph_deg = [0, 0, 90, -70];
th_deg = [0,90,180,270];
pol_deg = [0,0,0,0];

plots = 0; %boolean to turn plots on or off


%%%%%%%%%%%%% basis rotation %%%%%%%%%%%
th_deg = th_deg-45;
disp('Starting Calculation and Timer')
tic
disp('%%%%%%%%%%%%%%%Computing Lattice Parameters%%%%%%%%%%%%%%%')
[waveAmplitudes,deltaKsUnitless,deltaPhis,maxAmp]=GeneralLatticeComplexRepWithComponents(A,ph_deg,th_deg,pol_deg,plots);

%Make lattice manually...
if (0)
    disp('%%%%%%%%%%%%%%%%%% Manual Lattice Parameter Override %%%%%%%%%%%%')
    %Since these are made negative later
    waveAmplitudes = [0.5,0.5];
    deltaKsUnitless = [1 1;1 -1];
    deltaPhis = [0 0];
    [X,Y] = meshgrid(0:0.1:2,0:0.1:2);
    totalFromCompSum = zeros(size(X));
    for ii = 1:length(waveAmplitudes)
        %note when using the unitless k vectors to plot in untiless x as a
        %function of the lattice light wavelength need a 2pi factor
        compWaveMatrix = waveAmplitudes(ii)*cos(deltaKsUnitless(ii,1).*X*2*pi + deltaKsUnitless(ii,2).*Y*2*pi - deltaPhis(ii));
        totalFromCompSum = totalFromCompSum + compWaveMatrix;
    end
    maxAmp = max(max(totalFromCompSum))-min(min(totalFromCompSum));
%     figure
%     surf(X,Y,totalFromCompSum);
end
toc
disp("%%%%%%%%%%%%%%% Constructing Hamiltonian %%%%%%%%%%%%%%%")
tic
%since these are effectively indices (in fourier space), we need these to
%be rounded to integers. They pretty much already are to MATLAB precision

%NOTE, there is a bit of weirdness going on here in the sense that really I
%should be rescaling the basis vectors when I rotate the axes by 45 
%degrees. However, right now the round function is actually taking care of
%that. But, for other lattice configurations this could be broken. A note
%for fixing later... (I made this note around 9/10/2020)
deltaKsUnitless = round(deltaKsUnitless);
potentialDepth = 30; %in Er!!
waveAmplitudes = waveAmplitudes.*(potentialDepth./maxAmp);
%% Find Fourier Coefficients
%Size of the plane wave basis. I.E. the
%fourier coefficients will be included up to +- max_m. This means that the
%matrix size along one dimension will be 2*max_m + 1
max_m = 5;
mLength = 2*max_m + 1;
Vcoeff = zeros(mLength,mLength);
for jj = 1:length(waveAmplitudes)
    xKcomp = deltaKsUnitless(jj,1);
    yKcomp = deltaKsUnitless(jj,2);
    %I am making the index xKcomp+(max_m+1)because x/y K components could
    %be negative. Therefore, I want to make sure that the 'center' of the
    %matrix corresponds to the (0,0) coefficient. This is also related to
    %my above comment, since if the basis is too small, then these
    %arguments could be negative and matlab will complain :(
    Vcoeff(xKcomp+(max_m+1),yKcomp+(max_m+1)) = Vcoeff(xKcomp+(max_m+1),yKcomp+(max_m+1)) + waveAmplitudes(jj).*(exp(-1i*deltaPhis(jj)))./2;
    Vcoeff(-xKcomp+(max_m+1),-yKcomp+(max_m+1)) = Vcoeff(-xKcomp+(max_m+1),-yKcomp+(max_m+1)) + waveAmplitudes(jj).*(exp(1i*deltaPhis(jj)))./2;
end
%I have added this here since I think that Peter's code outputs the
%intensity
Vcoeff = -Vcoeff;

%% Plot the Potential (Optional, mainly for troubleshooting...)
figure
tic
trans_pott = fourier_to_real_fft(Vcoeff,max_m,50);
toc
contourf(repmat(trans_pott,5,5),10);
% surf(trans_pott)
% keyboard;
%% Create Hamiltonian
%now this is where things start to get different from the other code, given
%that we are considering only certain values of the quasimomentum for
%computing the wannier states

L = 7; %how many real space lattice sites to use. MAKE AN ODD NUMBER
max_qm = floor(L/2); %to make sure that we are always in the first BZ
qsize = 2*max_qm + 1;
zone_number = 1; %how many zones to plot (for extended zone picture)
[quasiX,quasiY] = meshgrid((-max_qm:max_qm)./L,(-max_qm:max_qm)./L);
quasiX = quasiX.';
quasiY = quasiY.';
%hammy is the hamiltonian. Needs to accommodate all of the information in
%components but in 2-d
hammy = zeros(mLength^2,mLength^2,qsize,qsize);
for ii = 1:mLength
    for jj = 1:mLength
        for kk = 1:mLength
            for ll = 1:mLength
                %to transform to the code, each of the indices actual
                %values in terms of the matlab indices are:
                % ii(actual) = ii - (max_m+1). This way the center of the
                % matrix corresponds to momentum index (0,0).
                if (ii==kk && jj==ll)
%                     size(components(ii,jj,kk,ll,:,:))
%                     size(((quasiX+ii-(max_m+1)).*(quasiX+ii-(max_m+1))))
                    hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:)=reshape(hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:),qsize,qsize)+((quasiX+ii-(max_m+1)).*(quasiX+ii-(max_m+1)))+((quasiY+jj-(max_m+1)).*(quasiY+jj-(max_m+1)));
                end
                if (abs(ii-kk) <= max_m && abs(jj-ll) <= max_m)
                      hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:) = hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:)+ Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
                end
            end
        end
    end
end

num_bands = 2;
eigvals = zeros(num_bands,num_bands,qsize,qsize);
eigvecs = zeros((mLength^2),num_bands,qsize,qsize);
toc
disp("%%%%%%%%%%%%Diagonalizing Hamiltonian %%%%%%%%%%%%%%%")
tic
for ii = 1:qsize
    for jj = 1:qsize
        %diagonalize for each point in the BZ
        [eigvecs(:,:,ii,jj),eigvals(:,:,ii,jj)] = eigs(hammy(:,:,ii,jj),num_bands,'smallestreal');
    end
end


%% Plot the band structure
toc
if (0)
    disp("%%%%%%%%%%%%%%% Plotting Band Structure %%%%%%%%%%%%%%%%")
    tic
    % fig1 = figure;
    figure;
    hold all;


    for kk = 1:num_bands
        surf(quasiX,quasiY,reshape(eigvals(kk,kk,:,:),qsize,qsize),'edgecolor','interp');
        xlabel('x quasimomentum, [$k_l$]','interpreter','latex');
        ylabel('y quasimomentum, [$k_l$]','interpreter','latex');
        zlab = ['Energy, [$E_R$]; lattice depth = ' num2str(potentialDepth) ' [$E_R$]'];
        zlabel(zlab, 'interpreter','latex');
    end

    toc
end
%% De-compactify the eigenvectors:
%bands to be caluclated. Really there are only 2 non-trivial sublattice
%minima
bands = [1 2];
disp("%%%%%%%%%%%%%%% De-Compactify Eigenvectors %%%%%%%%%%%%%%%%")
tic
%this is the eigenvector output from the eig function. Remember that we
%'compactified' the dimension here, so we want to make it back into the
%natural matrix form. This is really just for ease of use...
weightsMatrix = zeros(mLength,mLength,length(bands),qsize,qsize);
for kk = 1:length(bands)
    for ll = 1:qsize
        for mm = 1:qsize
            weightsMatrix(:,:,kk,ll,mm) = reshape(eigvecs(:,kk,ll,mm),mLength,mLength).';
        end
    end
end
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

% keyboard;
%perhaps here we compute all of the possible interations between the
%different wannier functions specified here?
%which x positions do we want to consider? Note these aren't the actual
%value of the real space x operator, but rather the inices of the latice
%unit cell. In this case, I want to take the one in the center as well as
%the one next to it so that I can compute tunelling element between them
x_positions = [ceil(L/2) (ceil(L/2)+1)];
%and the same story for the y positions
y_positions = [ceil(L/2) (ceil(L/2)+1)];
%which wannier funcs do we want? For this lattice we have two, and you
%could choose either one. Make sure that they are valid indices
sub_unit_funcs = [1 2]; 
grouped_eigenvecs = zeros(length(D2),L,length(x_positions),length(sub_unit_funcs));
grouped_eigenvals = zeros(L,length(x_positions),length(sub_unit_funcs));
ind = kmeans(real(D2),L); 
%note that I am NOT happy with the following solution, given that it does
%not account for the edge case that there are unequal numbers of elements
%within some different clusters, but in the case that they are, and things
%are sorted, then we can just sort the indices to make them start at 1
ind = sort(ind);
%note the kmeans doesn't sort them according to the values of the
%datapoints, therefore I will have to rearrange the groups myself
for ii = 1:length(x_positions)
    %now within each position we need to group by sub-minima index wannier
    %function
    pos_eigvals = D2(ind==x_positions(ii));
    pos_eigvecs = V2(:,ind==x_positions(ii));
    %note the sort here. Since the eigenvalues are already sorted, so will
    %pos_eigvals, and therefore I want the sub_ind (which represents the
    %band label of the wannier function) to have a consistent
    %interpretation
    sub_ind = sort(kmeans(real(pos_eigvals),length(bands)));
    for jj = 1:length(sub_unit_funcs)
        grouped_eigenvecs(:,:,ii,jj) = pos_eigvecs(:,sub_ind==sub_unit_funcs(jj));
        grouped_eigenvals(:,ii,jj) = pos_eigvals(sub_ind==sub_unit_funcs(jj));
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

wannier_states_in_bloch_basis = zeros(length(D2),length(x_positions),length(y_positions),length(sub_unit_funcs));
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

%% Generate Real-Space Plots
toc
disp("%%%%%%%%%%%%%%%Generate Real Space Plots%%%%%%%%%%%%%%%%")
tic


x_pos_index = 1;
y_pos_index = 1;
wannier_func_index = 1;
wannier_func_realspace = bloch_to_real_fft(wannier_states_in_bloch_basis(:,x_pos_index,y_pos_index,wannier_func_index),weightsMatrix,max_m,L,bands,0);

    
if (1)
    %Then just take away the imaginary part:
    wannier_func_realspace = wannier_func_realspace .* exp(-1i*angle(wannier_func_realspace(53,55))).*(-1);
    fontsize = 26;
    figure
    surf(-1.*real(wannier_func_realspace),'edgealpha',0.6);
    xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    ylabel('Y Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    title(['MLWF, $V_0 = $' num2str(potentialDepth) '$E_{R}$' ],'interpreter','latex','fontsize',fontsize);
    zlab = ['Re(Wannier Func)'];
    zlabel(zlab, 'interpreter','latex','fontsize',fontsize);
    figure
    surf(imag(wannier_func_realspace));
    xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    ylabel('Y Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    zlab = ['Im(Wannier Func)'];
    zlabel(zlab, 'interpreter','latex','fontsize',fontsize);
    figure
    surf(abs(wannier_func_realspace),'edgealpha',0.3);
%     surf(X,Y,abs(wannier_func_realspace));
    xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    ylabel('Y Pos., [$\lambda_l$]','interpreter','latex','fontsize',fontsize);
    zlab = ['Abs(Wannier)'];
    title('Wannier Func 1','interpreter','latex','fontsize',fontsize);
    zlabel(zlab, 'interpreter','latex','fontsize',fontsize);
    figure
    surf(angle(wannier_func_realspace));
    xlabel('X Pos., [$\lambda_l$]','interpreter','latex');
    ylabel('Y Pos., [$\lambda_l$]','interpreter','latex');
    zlab = ['Arg(Wannier Func)'];
    zlabel(zlab, 'interpreter','latex');
end

%% Calculate Bose Hubbard Parameters
toc
disp('%%%%%%%%%%%%%%%% Calculate Bose Hubbard Parameters %%%%%%%%%%%%%%%%%%%')
%note that now we need to re-construct the hamiltonian in the sub-space of
%the bloch states that we are trying to consider. Since we already found
%the eigenbasis of the hamiltonian, this is just a question of construcing
%a diagonal matrix with the eigenvalues of the hamiltonian already found

%note that we only need the collection of bands that are specified in the
%calculation earlier.

sub_space_hamiltonian = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
vectorized_eigvals = zeros(qsize.^2,length(bands));
for ii = bands
   vectorized_eigvals(:,ii) = reshape(eigvals(ii,ii,:,:),qsize.^2,1);
end
sub_space_hamiltonian(:,:) = diag(reshape(vectorized_eigvals,(qsize.^2)*length(bands),1));


%remove (as much as possible) the imaginary phase)
for ii = 1:length(x_positions)
    wannier_func_realspace = bloch_to_real_fft(wannier_states_in_bloch_basis(:,ii,y_pos_index,wannier_func_index),weightsMatrix,max_m,L,bands,0);
    wannier_states_in_bloch_basis(:,ii,y_pos_index,wannier_func_index) = wannier_states_in_bloch_basis(:,ii,y_pos_index,wannier_func_index).*exp(-1i*angle(wannier_func_realspace(53,55))).*(-1);
end
%calculate the J matrix element with the inner product of the hamiltonian
%and the wannier states in the bloch basis...
%in this case let's just calculate going one across in the x direction for
%the first wannier state
J = (wannier_states_in_bloch_basis(:,1,1,1)')*sub_space_hamiltonian*wannier_states_in_bloch_basis(:,2,1,1)


end



