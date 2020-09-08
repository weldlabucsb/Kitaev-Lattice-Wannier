function [] =  numerics_testing()
A = [1,1,0.6,0.5];

ph_deg = [0, 0, 90, -70];

th_deg = [0,90,180,270];
pol_deg = [0,0,0,0];
plots = 0; 
th_deg = th_deg-45;
[waveAmplitudes,deltaKsUnitless,deltaPhis,maxAmp]=GeneralLatticeComplexRepWithComponents(A,ph_deg,th_deg,pol_deg,plots);
disp("%%%%%%%%%%%%%%% Constructing Hamiltonian %%%%%%%%%%%%%%%")
tic
%since these are effectively indices (in fourier space), we need these to
%be rounded to integers. They pretty much already are to MATLAB precision
deltaKsUnitless = round(deltaKsUnitless);
potentialDepth = 40; %in Er!!
waveAmplitudes = waveAmplitudes.*(potentialDepth./maxAmp);

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

% keyboard;
%% Create Hamiltonian

L = 5; %this is the number of lattice sites along one direction to consider. The total number of sites to consider is L^2
max_qm = floor(L/2); %to make sure that we are always in the first BZ
qsize = 2*max_qm + 1;
zone_number = 1; %how many zones to plot (for extended zone picture)
[quasiX,quasiY] = meshgrid((-max_qm:max_qm)./L,(-max_qm:max_qm)./L);
%hammy is the hamiltonian. Needs to accommodate all of the information in
%components but in 2-d
hammy = zeros(mLength^2,mLength^2,qsize,qsize);
for ii = 1:mLength
    for jj = 1:mLength
        for kk = 1:mLength
            for ll = 1:mLength
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

num_bands = 2;
eigvals = zeros(num_bands,num_bands,qsize,qsize);
eigvecs = zeros((mLength^2),num_bands,qsize,qsize);
toc
disp("%%%%%%%%%%%%Diagonalizing Hamiltonian %%%%%%%%%%%%%%%")
tic
for ii = 1:qsize
    for jj = 1:qsize
        [eigvecs(:,:,ii,jj),eigvals(:,:,ii,jj)] = eigs(hammy(:,:,ii,jj),num_bands,'smallestreal');
    end
end
bands = [1 2];
%% De-compactify the eigenvectors:

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
            weightsMatrix(:,:,kk,ll,mm) = reshape(eigvecs(:,kk,ll,mm),mLength,mLength).';
        end
    end
end
%great! now we have a more interpretable to get the information that we
%need
%% Construct band projected position operators the analytic way...
toc 
disp("%%%%%%%%%%%% Construct Band-Projected Position Operators %%%%%%%%%%%%%%%%")
disp('Parallel operation only');
tic
R1 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
R2 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
parfor oo = 1:numel(R1)
    [row,col] = ind2sub([qsize.^2*length(bands),qsize.^2*length(bands)],oo);
    if (row >= col)
        [qy1,qx1,band1] = ind2sub([qsize,qsize,length(bands)],row);
        [qy2,qx2,band2] = ind2sub([qsize,qsize,length(bands)],col);
        [R1(oo),R2(oo)] = comp_elem_dingue(weightsMatrix,L,max_m,max_qm,band1,band2,qx1,qx2,qy1,qy2);
    end
end
R1 = R1 + (R1-diag(diag(R1)))';
R2 = R2 + (R2-diag(diag(R2)))';




toc 
disp("%%%%%%%%%%%% Construct Band-Projected Position Operators %%%%%%%%%%%%%%%%")
disp('Now with tensorized computation AND parallelization');
tic
R5 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
R6 = zeros(qsize.^2*length(bands),qsize.^2*length(bands));
parfor oo = 1:numel(R5)
    [row,col] = ind2sub([qsize.^2*length(bands),qsize.^2*length(bands)],oo);
    if (row >= col)
        [qy1,qx1,band1] = ind2sub([qsize,qsize,length(bands)],row);
        [qy2,qx2,band2] = ind2sub([qsize,qsize,length(bands)],col);
        [R5(oo),R6(oo)] = comp_elem_ok(weightsMatrix,L,max_m,max_qm,band1,band2,qx1,qx2,qy1,qy2);
    end
end
R5 = R5 + (R5-diag(diag(R5)))';
R6 = R6 + (R6-diag(diag(R6)))';
toc
keyboard;
end


function [x_element,y_element]=comp_elem_dingue(weightsMatrix,L,max_m,max_qm,ii,jj,kk,ll,mm,nn)
    x_element = 0;
    y_element = 0;
    m0 = ceil((L-1)./2)+0.5;
    mLength = 2*max_m + 1;
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
                            %note that I think there is a problem here with
                            %the 
                            x_element = x_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                .*((L^3)./2 - m0.*L^2);
                            y_element = y_element + conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn).*exp(-2*pi*1i*m0*(c1+c2))...
                                .*((L^3)./2 - m0.*L^2);
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
end

function [x_element,y_element]=comp_elem_ok(weightsMatrix,L,max_m,max_qm,ii,jj,kk,ll,mm,nn)
    m0 = ceil((L-1)./2)+0.5;
    mLength = 2*max_m + 1;
    braweightsMatrix = weightsMatrix(:,:,ii,kk,mm);
    ketweightsMatrix = weightsMatrix(:,:,jj,ll,nn);
    fourdweights = repmat(conj(braweightsMatrix),1,1,mLength,mLength);
    for pp = 1:mLength
        for rr = 1:mLength
            fourdweights(:,:,pp,rr) = fourdweights(:,:,pp,rr).*ketweightsMatrix(pp,rr);
        end
    end
    [oo,qq,pp,rr] = ndgrid(1:mLength,1:mLength,1:mLength,1:mLength);

%trying a quick re-definition since things were a little differently
%defined up there. 
%     n1prime = oo-max_m-1;
%     n2prime = qq-max_m-1;
%     n1 = pp-max_m-1;
%     n2 = rr-max_m-1;
%     m1prime = kk-max_qm-1;
%     m1 = ll-max_qm-1;
%     m2prime = mm-max_qm-1;
%     m2 = nn-max_qm-1;
    c1 = (oo-max_m-1-pp-max_m-1)+(kk-max_qm-1-ll-max_qm-1)./L;
    c2 = (qq-max_m-1-rr-max_m-1)+(mm-max_qm-1-nn-max_qm-1)./L;
    
    %first compute the c1 short integral. Use is nan to rectify the
    %removable singularities.
    c1short = (exp(1i*L*2*pi*c1)-1)./(1i*2*pi*c1);
    c1short(isnan(c1short)) = L;
    
    %next short c2 integral.
    c2short = (exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2);
    c2short(isnan(c2short)) = L;
    
    %next long c1 integral
    c1long = (exp(1i*L*2*pi*c1).*(1i*L*2*pi*c1-1)+1)./(-4*(pi^2)*(c1.*c1));
    c1long(isnan(c1long)) = L^2/2;
    
    %finally long c2 integral
    c2long = (exp(1i*L*2*pi*c2).*(1i*L*2*pi*c2-1)+1)./(-4*(pi^2)*(c2.*c2));
    c2long(isnan(c2long)) = L^2/2;
    
    %Great, now let's combine everything!
    
    x_element = (2*pi./L^2)*fourdweights.*exp(-1i*2*pi*m0*(c1+c2)).*(c1long.*c2short-m0.*c1short.*c2short);
    x_element = sum(x_element,'all');
    y_element = (2*pi./L^2)*fourdweights.*exp(-1i*2*pi*m0*(c1+c2)).*(c2long.*c1short-m0.*c1short.*c2short);
    y_element = sum(y_element,'all');
end

function [component_part_x,component_part_y,integral_part_x,integral_part_y,x_element,y_element]=comp_elem_dingue_testing(weightsMatrix,L,max_m,max_qm,ii,jj,kk,ll,mm,nn)
    m0 = ceil((L-1)./2)+0.5;
    mLength = 2*max_m + 1;
    x_element = 0;
    y_element = 0;
    component_part_x = zeros(mLength,mLength,mLength,mLength);
    component_part_y = zeros(mLength,mLength,mLength,mLength);
    integral_part_x = zeros(mLength,mLength,mLength,mLength);
    integral_part_y = zeros(mLength,mLength,mLength,mLength);
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
                    component_part = conj(weightsMatrix(oo,qq,ii,kk,mm)).*weightsMatrix(pp,rr,jj,ll,nn);
                    component_part_x(oo,pp,qq,rr) = component_part;
                    component_part_y(oo,pp,qq,rr) = component_part;
                    if c1 == 0
                        if c2 == 0
                            %note!!! I changed from the original code here
                            %to check
                            x_element = x_element + component_part.*exp(-2*pi*1i*m0*(c1+c2))...
                                .*((L^3)./2 - m0.*L^2);
                            integral_part_x(oo,pp,qq,rr) = exp(-2*pi*1i*m0*(c1+c2))...
                                .*((L^3)./2 - m0.*L^2);
                            y_element = y_element + component_part.*exp(-2*pi*1i*m0*(c1+c2))...
                                .*((L^3)./2 - m0.*L^2);
                            integral_part_y(oo,pp,qq,rr) = exp(-2*pi*1i*m0*(c1+c2))...
                                .*((L^3)./2 - m0.*L^2);
                            
                        else
                            x_element = x_element + component_part.*exp(-2*pi*1i*m0*(c1+c2))...
                                .*((exp(1i*L*2*pi*c2)-1)*((L^2)./2)./(1i*2*pi*c2)-m0*(exp(1i*L*2*pi*c2)-1)*L./(1i*2*pi*c2));
                            integral_part_x(oo,pp,qq,rr) = exp(-2*pi*1i*m0*(c1+c2))...
                                .*((exp(1i*L*2*pi*c2)-1)*((L^2)./2)./(1i*2*pi*c2)-m0*(exp(1i*L*2*pi*c2)-1)*L./(1i*2*pi*c2));
                            y_element = y_element + component_part.*exp(-2*pi*1i*m0*(c1+c2))...
                                .*(((exp(1i*L*2*pi*c2)*(1i*L*2*pi*c2-1)+1)./(-(2*pi*c2)^2))*L-m0*(exp(1i*L*2*pi*c2)-1)*L./(1i*2*pi*c2));
                            integral_part_y(oo,pp,qq,rr) = exp(-2*pi*1i*m0*(c1+c2))...
                                .*(((exp(1i*L*2*pi*c2)*(1i*L*2*pi*c2-1)+1)./(-(2*pi*c2)^2))*L-m0*(exp(1i*L*2*pi*c2)-1)*L./(1i*2*pi*c2));
                        end
                    else
                        if c2 == 0
                            x_element = x_element + component_part.*exp(-2*pi*1i*m0*(c1+c2))...
                                .*(((exp(1i*L*2*pi*c1)*(1i*L*2*pi*c1-1)+1)./(-(2*pi*c1)^2))*L-m0*(exp(1i*L*2*pi*c1)-1)*L./(1i*2*pi*c1));
                            integral_part_x(oo,pp,qq,rr) = exp(-2*pi*1i*m0*(c1+c2))...
                                .*(((exp(1i*L*2*pi*c1)*(1i*L*2*pi*c1-1)+1)./(-(2*pi*c1)^2))*L-m0*(exp(1i*L*2*pi*c1)-1)*L./(1i*2*pi*c1));
                            y_element = y_element + component_part.*exp(-2*pi*1i*m0*(c1+c2))...
                                .*((exp(1i*L*2*pi*c1)-1)*((L^2)./2)./(1i*2*pi*c1)-m0*(exp(1i*L*2*pi*c1)-1)*L./(1i*2*pi*c1));
                            integral_part_y(oo,pp,qq,rr) = exp(-2*pi*1i*m0*(c1+c2))...
                                .*((exp(1i*L*2*pi*c1)-1)*((L^2)./2)./(1i*2*pi*c1)-m0*(exp(1i*L*2*pi*c1)-1)*L./(1i*2*pi*c1));
                        else
                            x_element = x_element + component_part.*exp(-2*pi*1i*m0*(c1+c2))...
                                .*(((exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2)).*((exp(1i*L*2*pi*c1)*(1i*L*2*pi*c1-1)+1)./(-(2*pi*c1)^2))...
                                -m0*((exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2))*((exp(1i*L*2*pi*c1)-1)./(1i*2*pi*c1)));
                            integral_part_x(oo,pp,qq,rr) = exp(-2*pi*1i*m0*(c1+c2))...
                                .*(((exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2)).*((exp(1i*L*2*pi*c1)*(1i*L*2*pi*c1-1)+1)./(-(2*pi*c1)^2))...
                                -m0*((exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2))*((exp(1i*L*2*pi*c1)-1)./(1i*2*pi*c1)));
                            y_element = y_element + component_part.*exp(-2*pi*1i*m0*(c1+c2))...
                                .*(((exp(1i*L*2*pi*c1)-1)./(1i*2*pi*c1)).*((exp(1i*L*2*pi*c2)*(1i*L*2*pi*c2-1)+1)./(-(2*pi*c2)^2))...
                                -m0*((exp(1i*L*2*pi*c2)-1)./(1i*2*pi*c2))*((exp(1i*L*2*pi*c1)-1)./(1i*2*pi*c1)));
                            integral_part_y(oo,pp,qq,rr) = exp(-2*pi*1i*m0*(c1+c2))...
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
end