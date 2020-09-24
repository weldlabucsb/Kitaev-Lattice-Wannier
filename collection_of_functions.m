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