function [] =  numerics_testing()
A = [1,1,0.6,0.5];
ph_deg = [0, 0, 90, -70];
th_deg = [0,90,180,270];
pol_deg = [0,0,0,0];
plots = 0;
[waveAmplitudes,deltaKsUnitless,deltaPhis,maxAmp]=GeneralLatticeComplexRepWithComponents(A,ph_deg,th_deg,pol_deg,plots);
%% The Slow Way
%%%%%%%%%%%%%%%%%%%%%%%%%% The slow way %%%%%%%%%%%%%

max_m = 12;
qsize = 20;

tic 
deltaKsUnitless = round(deltaKsUnitless);
potentialDepth = 5; %in Er!!
waveAmplitudes = waveAmplitudes.*(potentialDepth./maxAmp);


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

%qusimomenta that you want
zone_number = 1; %how many extended zones to plot
[quasiX,quasiY] = meshgrid(linspace(-0.5*zone_number,0.5*zone_number,qsize),linspace(-0.5*zone_number,0.5*zone_number,qsize));
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
%                     components(ii,jj,kk,ll,:,:) = components(ii,jj,kk,ll,:,:) + Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
                      hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:) = hammy((mLength*(ii-1)+jj),(mLength*(kk-1)+ll),:,:)+ Vcoeff((ii-kk+(max_m+1)),(jj-ll+(max_m+1)));
                end
            end
        end
    end
end


eigvals = zeros((mLength^2),(mLength^2),qsize,qsize);
eigvecs = zeros((mLength^2),(mLength^2),qsize,qsize);
toc

disp("%%%%%%%%%%%%Diagonalizing Hamiltonian %%%%%%%%%%%%%%%")

for ii = 1:qsize
    for jj = 1:qsize
        [eigvecs(:,:,ii,jj),eigvals(:,:,ii,jj)] = eig(hammy(:,:,ii,jj));
    end
end

toc
keyboard;
%% Quick Way Now
%%%%%%%%%%%%%%%%%%%% the quick way %%%%%%%%%%%%%%%%%%%%
tic 
deltaKsUnitless = round(deltaKsUnitless);
potentialDepth = 5; %in Er!!
waveAmplitudes = waveAmplitudes.*(potentialDepth./maxAmp);

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
%qusimomenta that you want
zone_number = 1; %how many extended zones to plot
[quasiX,quasiY] = meshgrid(linspace(-0.5*zone_number,0.5*zone_number,qsize),linspace(-0.5*zone_number,0.5*zone_number,qsize));
%hammy is the hamiltonian. Needs to accommodate all of the information in
%components but in 2-d
% hammy = zeros(mLength^2,mLength^2,qsize,qsize);
hammy = cell(qsize);
hammy(:,:) = {sparse(mLength^2)};



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


num_bands = 5;


eigvals = zeros(num_bands,num_bands,qsize,qsize);
eigvecs = zeros((mLength^2),num_bands,qsize,qsize);
toc
disp("%%%%%%%%%%%%Diagonalizing Hamiltonian %%%%%%%%%%%%%%%")

for ii = 1:qsize
    for jj = 1:qsize
        [eigvecs(:,:,ii,jj),eigvals(:,:,ii,jj)] = eigs(hammy(:,:,ii,jj),num_bands,'smallestabs');
    end
end


toc
keyboard;


end