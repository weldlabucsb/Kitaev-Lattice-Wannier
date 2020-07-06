function[] = lattice_wannier()



%COMPUTEBANDSTRUC Compute 2-d band structure of optical lattice 
%   Using the potential output by the general lattice for N beams code, I
%   will use this program to compute the electronic band structure of the
%   system. Note, this needs to be in the same directory (or path) as the
%   general lattice complex function given that the latter needs to be
%   called to generate the potential (specifically the underlying
%   sinusoidal representation)

%% generate the lattice potential
%using the 6 cosine components of the lattice potential generated as output
%from Peter's code.
% A = [1,1,0,0];
A = [1,1,0.6,0.5];
% ph_deg = [0, 0, 0, 0];
ph_deg = [0, 0, 90, -70];
% th_deg = [0,90,180,270];
th_deg = [0,90,180,270];
pol_deg = [0,0,0,0];
% pol_deg = [0,0,0,0];
% A = [1,1,1,1];
% ph_deg = [0, 0,250,60];
% th_deg = [0,90,180,270];
% pol_deg = [0,0,0,0];
plots = 0; %boolean to turn plots on or off

%%%%%%%%%% Alternate lattice params here
% A = [1,1,1,1];
% ph_deg = [0, 0,250,60];
% th_deg = [0,90,180,270];
% pol_deg = [0,0,0,0];
% plots = 1;
%%%%%%%%%%%%

%%%%%%%%%%%%% basis rotation %%%%%%%%%%%
th_deg = th_deg-45;


[waveAmplitudes,deltaKsUnitless,deltaPhis,maxAmp]=GeneralLatticeComplexRepWithComponents(A,ph_deg,th_deg,pol_deg,plots);


%Make lattice manually...
if (0)
    %Since these are made negative later
    waveAmplitudes = [-1/2,-1/2];
    deltaKsUnitless = [1 1;1 -1];
    deltaPhis = [0 0];
    maxAmp = 2;
    [X,Y] = meshgrid(0:0.1:2,0:0.1:2);
    totalFromCompSum = zeros(size(X));
    for ii = 1:length(waveAmplitudes)
        %note when using the unitless k vectors to plot in untiless x as a
        %function of the lattice light wavelength need a 2pi factor
        compWaveMatrix = waveAmplitudes(ii)*cos(deltaKsUnitless(ii,1).*X*2*pi + deltaKsUnitless(ii,2).*Y*2*pi - deltaPhis(ii));
        totalFromCompSum = totalFromCompSum + compWaveMatrix;
    end
    maxAmp = max(max(totalFromCompSum))-min(min(totalFromCompSum));
    figure
    surf(X,Y,totalFromCompSum);
end

%since these are effectively indices (in fourier space), we need these to
%be rounded to integers. They pretty much already are to MATLAB precision
deltaKsUnitless = round(deltaKsUnitless);
potentialDepth = 10; %in Er!!
waveAmplitudes = waveAmplitudes.*(potentialDepth./maxAmp);
%% Find the Complex Exponential Coefficients
%Effectively I just want the coefficients of the complex fourier series

%BIG varible here. This is the size of the plane wave basis. I.E. the
%fourier coefficients will be included up to +- max_m. This means that the
%tensor size will be 2*(max_m)+1. Note that this MUST be large enough to
%include the highest frequency cosine that will be part of the potential.
%This is crucial since you could theoretically cut the matrix off at any
%point where you would be happy with the accuracy in finding the
%eigenvalues and vectors, but the KVectors could be arbitrarily large. In
%this case they are not (<=2), but in principle they could be quite large.
max_m = 4;
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

% keyboard;
%% Range of Quasimomenta to use
%now this is where things start to get different from the other code, given
%that we are considering only certain values of the quasimomentum for
%computing the wannier states

L = 5; %this is the number of lattice sites along one direction to consider. The total number of sites to consider is L^2
max_qm = floor(L/2); %to make sure that we are always in the first BZ
qsize = 2*max_qm + 1;
zone_number = 1; %how many zones to plot (for extended zone picture)
[quasiX,quasiY] = meshgrid((-max_qm:max_qm)./L,(-max_qm:max_qm)./L);
%hammy is the hamiltonian. Needs to accommodate all of the information in
%components but in 2-d
hammy = zeros(mLength^2,mLength^2,qsize,qsize);
tic
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
%                     components(ii,jj,kk,ll,:,:) = components(ii,jj,kk,ll,:,:) + Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
                      hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:) = hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:)+ Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
                end
            end
        end
    end
end
toc

eigvals = zeros((mLength^2),(mLength^2),qsize,qsize);
eigvecs = zeros((mLength^2),(mLength^2),qsize,qsize);

for ii = 1:qsize
    for jj = 1:qsize
        [eigvecs(:,:,ii,jj),eigvals(:,:,ii,jj)] = eig(hammy(:,:,ii,jj));
    end
end


%% Plot the band structure

% fig1 = figure;
figure;
hold all;

num_bands = 2;
for kk = 1:num_bands
    surf(quasiX,quasiY,reshape(eigvals(kk,kk,:,:),qsize,qsize),'edgecolor','interp');
    xlabel('x quasimomentum, [$k_l$]','interpreter','latex');
    ylabel('y quasimomentum, [$k_l$]','interpreter','latex');
    zlab = ['Energy, [$E_R$]; lattice depth = ' num2str(potentialDepth) ' [$E_R$]'];
    zlabel(zlab, 'interpreter','latex');
end

%bands to be caluclated. I think that we need 4 for the 4 sub-cell minima
bands = [1 2];

%% De-compactify the eigenvectors:
%this is the eigenvector output from the eig function. Remember that we
%'compactified' the dimension here, so we want to make it back into the
%natural matrix form. This is really just for ease of use...
weightsMatrix = zeros(mLength,mLength,length(bands),qsize,qsize);
for ii = 1:mLength
    for jj = 1:mLength
        for kk = 1:length(bands)
            for ll = 1:qsize
                for mm = 1:qsize
                    weightsColumn = eigvecs(:,kk,ll,mm);
                    weightsMatrix(ii,jj,kk,ll,mm) = weightsColumn((mLength*(ii-1)+jj));
                end
            end
        end
    end
end
%great! now we have a more interpretable to get the information that we
%need
%% Construct band projected position operators the analytic way...
%This is using the analytic expression given as a sum in the document. It
%is still quite long :/

m0 = ceil((L-1)./2)+0.5;

R1 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
R2 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
for ii = bands
    for jj = bands
        disp(['ii is ' num2str(ii) ', jj is ' num2str(jj)])
        for kk = 1:qsize
            for ll = 1:qsize
                for mm = 1:qsize
                    for nn = 1:qsize
                        x_element = 0;
                        y_element = 0;
                        
                        for oo = 1:mLength
                            for pp = 1:mLength
                                for qq = 1:mLength
                                    for rr = 1:mLength
                                        n1prime = oo-max_m-1;
                                        n1 = pp-max_m-1;
                                        n2prime = qq-max_m-1;
                                        n2 = rr-max_m-1;
                                        m1prime = kk-max_qm-1;
                                        m1 = ll-max_qm-1;
                                        m2prime = mm-max_qm-1;
                                        m2 = nn-max_qm-1;
                                        c1 = (n1prime-n1)+(m1prime-m1)./L;
                                        c2 = (n2prime-n2)+(m2prime-m2)./L;
                                        if c1 == 0
                                            if c2 == 0
                                                x_element = x_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                                    .*((L^3)./3 - m0.*L^2);
                                                y_element = y_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                                    .*((L^3)./3 - m0.*L^2);
                                            else
                                                x_element = x_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                                    .*((exp(1i*L*2*pi*c2)-1)*((L^2)./2)./(1i*2*pi*c2)-m0*(exp(1i*L*2*pi*c2)-1)*L./(1i*2*pi*c2));
                                                y_element = y_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                                    .*(((exp(1i*L*2*pi*c2)*(1i*L*2*pi*c2-1)+1)./(-(2*pi*c2)^2))*L-m0*(exp(1i*L*2*pi*c2)-1)*L./(1i*2*pi*c2));
                                            end
                                        else
                                            if c2 == 0
                                                x_element = x_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                                    .*(((exp(1i*L*2*pi*c1)*(1i*L*2*pi*c1-1)+1)./(-(2*pi*c1)^2))*L-m0*(exp(1i*L*2*pi*c1)-1)*L./(1i*2*pi*c1));
                                                y_element = y_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                                    .*((exp(1i*L*2*pi*c1)-1)*((L^2)./2)./(1i*2*pi*c1)-m0*(exp(1i*L*2*pi*c1)-1)*L./(1i*2*pi*c1));
                                            else
                                                x_element = x_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                                    .*(((exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2)).*((exp(1i*L*2*pi*c1)*(1i*L*2*pi*c1-1)+1)./(-(2*pi*c1)^2))...
                                                    -m0*((exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2))*((exp(1i*L*2*pi*c1)-1)./(1i*2*pi*c1)));
                                                y_element = y_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                                    .*(((exp(1i*L*2*pi*c1)-1)./(1i*2*pi*c1)).*((exp(1i*L*2*pi*c2)*(1i*L*2*pi*c2-1)+1)./(-(2*pi*c2)^2))...
                                                    -m0*((exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2))*((exp(1i*L*2*pi*c1)-1)./(1i*2*pi*c1)));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                         x_element = x_element*(2*pi./L^2);
                         y_element = y_element*(2*pi./L^2);
                        R1((ii-1)*(qsize)^2 + (kk-1)*qsize + mm,(jj-1)*qsize^2+(ll-1)*qsize+nn) = x_element;
                        R2((ii-1)*(qsize)^2 + (kk-1)*qsize + mm,(jj-1)*qsize^2+(ll-1)*qsize+nn) = y_element;
%                         toc
                    end
                end
            end
        end
    end
    toc
end
toc

%% Diagonalize the Band Projected Position Operators.
%Technically only one needs to be fully diagonalized first, since after we
%do this we will only need to diagonalize the other one in a
%semi-degenerate subspace of the first one

% e2 = eig(R2);
% e1 = eig(R1);
[V2,D2] = eig(R2);
E2 = diag(D2);
figure;
grouped_eigenvecs = zeros(1,L);
grouped_eigenvals = zeros(length(V2(:,1)),L);
jj = 1;
for ii = 1:length(E2)
    if (E2(ii) > -6.1) && (E2(ii) < -5.9)
        grouped_eigenvecs(jj) = E2(ii);
        grouped_eigenvals(:,jj) = V2(:,ii);
        jj = jj+1;
    end
end
jj %just a double check that a group of 5 came out


% plot(sort(real(e1)));
% title('xoperator eigenvalues');
% figure;
% plot(sort(real(e2)));
% title('yoperator eigenvalues');
% keyboard

end

function [U] = bloch_func(weights,max_m,xmat,ymat)
%NOTE: we have to be a bit careful when using meshgrid here. Unfortunately
%the A(i,j) matrix indexing convention conflicts with our traditional idea
%of (x,y) pairs. For (i,j), i represents a vertical displacement, but for (x,y)
%the first index represents horizontal displacement. If this isn't clear,
%try testing a simple example of meshgrid, or check MATLAB documentation
%site. Long story short we can just transpose the weights matrix to take
%this into account
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

function [U] = bloch_wave(eigvecs,max_m,xmat,ymat,quasiX,quasiY,bands)
%Now the goal of this function is to actually return the eigenfunctions of
%the Hamiltonian. The other ones are the bloch functions (that have the
%periodicity of the lattice)

%This is the grid of higher momenta in general.
[mxMat,myMat] = meshgrid(-max_m:max_m,-max_m:max_m);
mLength = 2*max_m + 1;
%Note that the first two indices are the real space representation of the
%bloch waves. The next two indices are the quasimomenta. The final one is
%the band index. I think that we should only need the first four, since we
%have four minima in a single lattice cell
U = zeros(length(xmat(:,1)),length(ymat(:,1)),length(bands),length(quasiX(:,1)),length(quasiY(:,1)));
for kk = bands
    %please see the note in the bloch function function about why this is
    %weightsMatrix(jj,ii) as opposed to (ii,jj).
    for ll = 1:length(quasiX(:,1))
        for mm = 1:length(quasiY(:,1))
            
            weightsColumn  = eigvecs(:,kk,ll,mm);
            weightsMatrix = zeros(mLength,mLength);
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