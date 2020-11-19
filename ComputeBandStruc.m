function [] = ComputeBandStruc()
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
plots = 1; %boolean to turn plots on or off

% A = [1,1,1,1];
% ph_deg = [0, 0,250,60];
% th_deg = [0,90,180,270];
% pol_deg = [0,0,0,0];
% plots = 1;


[waveAmplitudes,deltaKsUnitless,deltaPhis,maxAmp]=GeneralLatticeComplexRepWithComponents(A,ph_deg,th_deg,pol_deg,plots);
%square lattice params
if (1)
    waveAmplitudes = [1,1];
    deltaKsUnitless = [1 0;0 1];
    deltaPhis = [0 0];
    maxAmp = 2;
    [X,Y] = meshgrid(-4:0.1:4,-4:0.1:4);
    totalFromCompSum = zeros(size(X));
    for ii = 1:length(waveAmplitudes)
        compWaveMatrix = waveAmplitudes(ii)*cos(deltaKsUnitless(ii,1).*X + deltaKsUnitless(ii,2).*Y - deltaPhis(ii));
        totalFromCompSum = totalFromCompSum + compWaveMatrix;
    end
    figure
    surf(X,Y,totalFromCompSum);
end
%since these are effectively indices (in fourier space), we need these to
%be rounded to integers. They pretty much already are to MATLAB precision
deltaKsUnitless = round(deltaKsUnitless);
potentialDepth = 4; %in Er!!
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
max_m = 9;
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
%% Explicitly Create the Hamiltonian as a Matrix
% this may look simple yet confusing. There is a bit going on here.
% Hopefully the document can help explain
qsize = 15;
%qusimomenta that you want
zone_number = 1; %how many extended zones to plot
[quasiX,quasiY] = meshgrid(linspace(-0.5*zone_number,0.5*zone_number,qsize),linspace(-0.5*zone_number,0.5*zone_number,qsize));
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
% for ii = 1:mLength
%     for jj = 1:mLength
%         for kk = 1:mLength
%             for ll = 1:mLength
%                 hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:) = components(ii,jj,kk,ll,:,:);
%             end
%         end
%     end
% end
eigvals = zeros((mLength^2),(mLength^2),qsize,qsize);
eigvecs = zeros((mLength^2),(mLength^2),qsize,qsize);
for ii = 1:qsize
    for jj = 1:qsize
        [eigvecs(:,:,ii,jj),eigvals(:,:,ii,jj)] = eig(hammy(:,:,ii,jj));
    end
end
toc
%% Plot the band structure
% keyboard;
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

keyboard
% keyboard
%% Plot the bloch functions
%The bloch functions are themselves a function of the quasimomentum, so
%select a mometum point here to look at. 
qxIndex = 3;
qyIndex = qxIndex;
bandNumber = 2;
%this is the eigenvector output from the eig function. Remember that we
%'compactified' the dimension here, so we want to make it back into the
%natural matrix form
weightsColumn  = eigvecs(:,bandNumber,qxIndex,qyIndex);
weightsMatrix = zeros(mLength,mLength);
for ii = 1:mLength
    for jj = 1:mLength
        weightsMatrix(ii,jj) = weightsColumn((mLength*(ii-1)+jj));
    end
end
maxval = 3;
xMin = 0;
xMax = maxval;  %units of lambda
yMin = 0;
yMax = maxval;  
numpoints = 70;
[X,Y] = meshgrid(linspace(xMin,xMax,numpoints),linspace(yMin,yMax,numpoints));
U = bloch_func(weightsMatrix,max_m,X,Y);
toc
figure
surf(X,Y,abs(U));
xlabel('X Pos., [$\lambda_l$]','interpreter','latex');
ylabel('y Pos., [$\lambda_l$]','interpreter','latex');
zlab = ['Bloch Fxn. [$|u_{n = ' num2str(bandNumber)' '}|, \mathbf{q} = (' num2str(quasiX(qxIndex,qyIndex)) ',' num2str(quasiY(qxIndex,qyIndex)) ')$'];
zlabel(zlab, 'interpreter','latex');
figure
surf(X,Y,real(U));
xlabel('X Pos., [$\lambda_l$]','interpreter','latex');
ylabel('y Pos., [$\lambda_l$]','interpreter','latex');
zlab = ['Bloch Fxn. [$RE(u_{n = ' num2str(bandNumber)' '}), \mathbf{q} = (' num2str(quasiX(qxIndex,qyIndex)) ',' num2str(quasiY(qxIndex,qyIndex)) ')$'];
zlabel(zlab, 'interpreter','latex');


%% Plot Bloch Waves (I.E. The actual eigenvalues)
bands = [1 2 3 4 5];
U = bloch_wave(eigvecs,max_m,X,Y,quasiX,quasiY,bands);
figure
surf(X,Y,abs(U(:,:,bandNumber,qxIndex,qyIndex)));
xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',30);
ylabel('y Pos., [$\lambda_l$]','interpreter','latex','fontsize',30);
zlab = ['Bloch Wave [$|u_{n = ' num2str(bandNumber)' '}|, \mathbf{q} = (' num2str(quasiX(qxIndex,qyIndex)) ',' num2str(quasiY(qxIndex,qyIndex)) ')$'];
zlabel(zlab, 'interpreter','latex','fontsize',30);

figure
surf(X,Y,real(U(:,:,bandNumber,qxIndex,qyIndex)));
xlabel('X Pos., [$\lambda_l$]','interpreter','latex');
ylabel('y Pos., [$\lambda_l$]','interpreter','latex');
zlab = ['Bloch Wave [$RE(u_{n = ' num2str(bandNumber)' '}, \mathbf{q} = (' num2str(quasiX(qxIndex,qyIndex)) ',' num2str(quasiY(qxIndex,qyIndex)) ')$'];
zlabel(zlab, 'interpreter','latex');
% keyboard;

figure
surf(X,Y,abs(U(:,:,bandNumber,qxIndex,qyIndex)));
xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',30);
ylabel('y Pos., [$\lambda_l$]','interpreter','latex','fontsize',30);
zlab = ['$|\Psi(\mathbf(r))|$'];
zlabel(zlab, 'interpreter','latex','fontsize',30);
title('Delocalized State','interpreter','latex','fontsize',30);

figure
surf(X,Y,abs(U(:,:,1,qxIndex,qyIndex)));
xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',30);
ylabel('y Pos., [$\lambda_l$]','interpreter','latex','fontsize',30);
zlab = ['$|\Psi(\mathbf(r))|$'];
zlabel(zlab, 'interpreter','latex','fontsize',30);
title('Delocalized State','interpreter','latex','fontsize',30);

figure
surf(X,Y,abs(U(:,:,4,qxIndex,qyIndex)));
xlabel('X Pos., [$\lambda_l$]','interpreter','latex','fontsize',30);
ylabel('y Pos., [$\lambda_l$]','interpreter','latex','fontsize',30);
zlab = ['$|\Psi(\mathbf(r))|$'];
zlabel(zlab, 'interpreter','latex','fontsize',30);
title('Delocalized State','interpreter','latex','fontsize',30);


%% Construct band projected position operators
%The bloch wave function is designed to be more general and make bloch
%waves for all quasimomenta values for the given list of specified bands

%bands to be caluclated. I think that we need 4 for the 4 sub-cell minima
% bands = [1 2 3];
% toc
% U = bloch_wave(eigvecs,max_m,X,Y,quasiX,quasiY,bands);
% toc
% 
% %First let us make the projected x position operator:
% R1 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
% R2 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));   
% for ii = bands
%     for jj = bands
%         for kk = 1:qsize
%             for ll = 1:qsize
%                 for mm = 1:qsize
%                     for nn = 1:qsize
%                         x_part = conj(U(:,:,ii,kk,mm)).*X.*U(:,:,jj,ll,nn);
%                         y_part = conj(U(:,:,ii,kk,mm)).*Y.*U(:,:,jj,ll,nn);
%                         x_element = sum(x_part,'all')*(xMax-xMin)*(yMax-yMin)./(numpoints.^2);
%                         y_element = sum(y_part,'all')*(xMax-xMin)*(yMax-yMin)./(numpoints.^2);
%                         R1((ii-1)*(qsize)^2 + (kk-1)*qsize + mm,(jj-1)*qsize^2+(ll-1)*qsize+nn) = x_element;
%                         R2((ii-1)*(qsize)^2 + (kk-1)*qsize + mm,(jj-1)*qsize^2+(ll-1)*qsize+nn) = y_element;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% toc
% e1 = eig(R1);
% e2 = eig(R2);
% figure()
% plot(real(e1))
% figure
% plot(real(e2))
% keyboard
% % Just a quick normalization sanity check. MATLAB should automatically
% normalize the output eigenvectors. 
% % disp('now doing integral');
% % probdens = abs(U).*abs(U);
% % integral = 0;
% % for ii = 1:length(U(1,:))
% %     for jj = 1:length(U(:,1))
% %         integral = integral + (xMax./numpoints)^2*(probdens(ii,jj));
% %     end
% % end
% % disp('integral is')
% % integral

%% Now the analytic method of finding the wannier functions
%There is a lot of detail about this in the theory document that I have
%typed up. Feel free to give that one a read...




%% Try to Construct the Wannier States
%I certainly think that I have some thinking to do here about the best way
%to construct these states. Reading the thesis we need to adjust the phases
%to get maximum localization. Need to do some thinking about this
toc

%wannier states calculated separately for each band



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
% keyboard;
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

