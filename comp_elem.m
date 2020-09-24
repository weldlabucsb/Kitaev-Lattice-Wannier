function [x_element,y_element]=comp_elem(weightsMatrix,L,max_m,max_qm,ii,jj,kk,ll,mm,nn)
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
    c1 = (oo-pp)+(kk-ll)./L;
    c2 = (qq-rr)+(mm-nn)./L;
    
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
    
    %Now to combine everything, and sum:
    x_element = (2*pi./L^2)*fourdweights.*exp(-1i*2*pi*m0*(c1+c2)).*(c1long.*c2short-m0.*c1short.*c2short);
    x_element = sum(x_element,'all');
    y_element = (2*pi./L^2)*fourdweights.*exp(-1i*2*pi*m0*(c1+c2)).*(c2long.*c1short-m0.*c1short.*c2short);
    y_element = sum(y_element,'all');

end