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
A = [1,1,0.6,0.5];
ph_deg = [0, 0, 90, -70];
th_deg = [0,90,180,270];
pol_deg = [0,0,0,0];
plots = 0; %boolean to turn plots on or off
[waveAmplitudes,deltaKsUnitless,deltaPhis,maxAmp]=GeneralLatticeComplexRepWithComponents(A,ph_deg,th_deg,pol_deg,plots);
%since these are effectively indices (in fourier space), we need these to
%be rounded to integers. They pretty much already are to MATLAB precision
deltaKsUnitless = round(deltaKsUnitless);
potentialDepth = 3; %in Er!!
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
max_m = 5;
mLength = 2*max_m + 1;
Vcoeff = zeros(mLength,mLength);
for jj = 1:6
    xKcomp = deltaKsUnitless(jj,1);
    yKcomp = deltaKsUnitless(jj,2);
    %I am making the index xKcomp+(max_m+1)because x/y K components could
    %be negative. Therefore, I want to make sure that the 'center' of the
    %matrix corresponds to the (0,0) coefficient. This is also related to
    %my above comment, since if the basis is too small, then these
    %arguments could be negative and matlab will complain :(
    Vcoeff(xKcomp+(max_m+1),yKcomp+(max_m+1)) = Vcoeff(xKcomp+(max_m+1),yKcomp+(max_m+1)) + waveAmplitudes(jj)*(exp(-1i*deltaPhis(jj)))./2;
    Vcoeff(-xKcomp+(max_m+1),-yKcomp+(max_m+1)) = Vcoeff(xKcomp+(max_m+1),yKcomp+(max_m+1)) + waveAmplitudes(jj)*(exp(1i*deltaPhis(jj)))./2;
end
%% Explicitly Create the Hamiltonian as a Matrix
% this may look simple yet confusing. There is a bit going on here.
% Hopefully the document can help explain

qsize = 10;
%qusimomenta that you want
[quasiX,quasiY] = meshgrid(linspace(-1,1,qsize),linspace(-1,1,qsize));
size(quasiX)
components = zeros(mLength,mLength,mLength,mLength,qsize,qsize);
%hammy is the hamiltonian. Needs to accommodate all of the information in
%components but in 2-d
hammy = zeros(mLength^2,mLength^2,qsize,qsize);
for ii = 1:mLength
    for jj = 1:mLength
        for kk = 1:mLength
            for ll = 1:mLength
                %to transform to the code, each of the indices actual
                %values in terms of the matlab indices are:
                % ii(actual) = ii - (max_m+1)
                if (ii==kk && jj==ll)
%                     size(components(ii,jj,kk,ll,:,:))
%                     size(((quasiX+ii-(max_m+1)).*(quasiX+ii-(max_m+1))))
                    components(ii,jj,kk,ll,:,:)=((quasiX+ii-(max_m+1)).*(quasiX+ii-(max_m+1)))+((quasiY+jj-(max_m+1)).*(quasiY+jj-(max_m+1)));
                end
                if (abs(ii-kk) <= max_m && abs(jj-ll) <= max_m)
%                     components(ii,jj,kk,ll,:,:) = components(ii,jj,kk,ll,:,:) + Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
                      components(ii,jj,kk,ll,:,:) = Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
                end
            end
        end
    end
end


%% Casting in 2-d Matrix Form
%now is the time we need to make this a 2-dimensional matrix to find the
%eigenvalues. I'm not going to do this in the most intuitive way, just the
%quickest. 
for ii = 1:mLength
    for jj = 1:mLength
        for kk = 1:mLength
            for ll = 1:mLength
                hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:) = components(ii,jj,kk,ll,:,:);
            end
        end
    end
end

eigs = zeros((mLength^2),qsize,qsize);
for ii = 1:qsize
    for jj = 1:qsize
        e = eig(hammy(:,:,ii,jj));
        eigs(:,ii,jj) = e;
    end
end

%% Plot the band structure

% fig1 = figure;
hold all;

num_bands = 1;
for kk = 1:num_bands
    surf(quasiX,quasiY,reshape(real(eigs(kk,:,:)),10,10));
end
% plotmat = zeros(num_bands,qsize,qsize);
% for kk = 1:num_bands
%     for ii = 1:qsize
%         for jj = 1:qsize
%             plotmat(kk,ii,jj) = real(eigs(kk,ii,jj));
%         end
%     end
% end
% plotting_matrix = zeros(qsize,qsize);
% for kk = 1:num_bands
%     plotting_matrix = reshape(plotmat(kk,:,:),10,10);
%     surf(quasiX,quasiY,plotting_matrix);
% end

end


