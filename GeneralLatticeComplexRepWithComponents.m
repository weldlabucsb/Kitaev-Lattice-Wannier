function [waveAmplitudes,deltaKsUnitless,deltaPhis]=GeneralLatticeComplexRepWithComponents()
    %% THIS FUNCTION IS USED TO MODEL A GENERAL 2D OPTICAL POTENTIAL
        % Presently it is coded for multipled interfering laser beams (plane waves)
        % Enables setting of wave amplitude, wave phase, wave propagation
        % angle, and wave polarization angle (in vs out of phase)'
        %
        % CARE MUST BE TAKEN IF ONE WANTS TO EXTEND TO
        % MULTIPLE WAVES OF DIFFERENT WAVELENGTH.  
        %
        % Extention possible to waves with elliptical polarization,
        % but not implemented here.
        
    close all

    %% Initialize Position Space Mesh and Time Mesh
    %units of lamba, the lattice light wavelength
    xMin = 0;
    xMax = 1;  %units of lambda
    yMin = 0;
    yMax = 1;   %units of lambda
    xStep = 0.01;
    yStep = 0.01;
    
    xVals = xMin:xStep:xMax;
    yVals = yMin:yStep:yMax;

    [X,Y] = meshgrid(xVals,yVals);
    
    kMag = 2*pi;  %Magnitude of kvectors.



    %% Set parameters for the laser beams

    
    % Electric field amplitudes
    A = [1,1,0.6,0.5];

    
    % Light beam phases
    ph_deg = [0, 0, 90, -70];
    ph = (pi/180)*ph_deg;

    
    %Planewave propagation direction angles (0 = positive x direction, 90deg = positive y direction)
    th_deg = [0,...
        90,...
        180,...
        270];

    th = (pi/180)*th_deg;


    % Polarization angles of the beams (0 = out of plane;  90deg = in plane)
    pol_deg = [0,0,0,0];
    pol = (pi/180)*pol_deg;


    %% Calculations of interference 
    numWaves = max([length(A),length(ph),length(th_deg),length(pol)]);

    E = zeros(size(X,1),size(X,2),numWaves); %initialize efield:  E(x,y,waveIndex)

    % Calculate and save efield complex representation as a function of space 
    % for each beam separately
    
    EzAmp = zeros(1,numWaves);
    Epara = zeros(1,numWaves);
    for qq = 1:numWaves
        kDotR = kMag*(cos(th(qq))*X + sin(th(qq))*Y); % Matrix of kDotR evalutated at all points in grid
        E(:,:,qq) = A(qq)*exp(1i*(kDotR -ph(qq)));  % Complex representation of eField Magnitude (i.e. not as a vector)
        EzAmp(qq) = A(qq)*cos(pol(qq));
        Epara(qq) = A(qq)*sin(pol(qq));
        
    end

    % Initialize components of the E field:   
    % complex vector E(x,y,waveIndex) = Ex*xhat + Ey*yhat +  Ez*zhat
    Ex = zeros(size(E));
    Ey = zeros(size(E));
    Ez = zeros(size(E));
    for qq = 1:numWaves
        Ex(:,:,qq) = E(:,:,qq)*(-sin(th(qq))*sin(pol(qq)));
        Ey(:,:,qq) = E(:,:,qq)*(cos(th(qq))*sin(pol(qq)));
        Ez(:,:,qq) = E(:,:,qq)*(cos(pol(qq)));
    end
    % Calculating the total x,y,z components of electric field by summing the
    % cartesian vector components of the complex electric field
    % representation
    
    % E_totx(x,y)
    E_totx = sum(Ex,3);  % sums along index 4 (the waveIndex)
    E_toty = sum(Ey,3);
    E_totz = sum(Ez,3);

    % Calculating the average intensity over one cycle.  Iavg is a function of all space.
    % This follows from an easy generalization of Zangwill's 
    %"A Time-Averaging Theorem" equation (1.139)
    % [Andrew Zangwill, MODERN ELECTRODYNAMICS, Cambridge University Press 2012, ISBN 978-0-521-89697-9]
    avgI = 0.5*(E_totx.*conj(E_totx) + E_toty.*conj(E_toty) + E_totz.*conj(E_totz));  


    %% Creating figure of average intensity as a function of 

    %Text to show paramters at the top of the figure
    EAmps = ['E-FieldAmps = [ ', num2str(A) , ' ]'];
    beamDirections = ['Beam Directions = [ ', num2str(th_deg) , ' ] (deg)'];
    beamPolarizationAngles = ['Beam Polarization Angles = [ ', num2str(pol_deg) , ' ] (degs)'];
    beamPhasesAngles = ['Beam Phases = [ ', num2str(ph_deg) , ' ] (deg)'];


    % Figure as a contour plot with N levels
    N=10;
    f1 = figure;
%     figHeight = 700;
%     set(f1,'Position',[20,950-figHeight,800,figHeight])
    set(f1,'color','w')
    contourf(X,Y,(avgI),N)
    xlabel('X Position, [$\lambda$]','interpreter','latex')
    ylabel('Y Position, [$\lambda$]','interpreter','latex')
    colorbar
    colormap(parula)
    set(f1,'Units','Pixels');

    % Adding parameters to the top of the figure
    uicontrol(f1,'style','text','Units','Pixels',...
        'Position',[10 f1.Position(4)-20 800 15],'String',[EAmps,';  ',beamPolarizationAngles,';  '], ...
        'HorizontalAlignment', 'Left','Background','w'); 
    uicontrol(f1,'style','text','Units','Pixels',...
        'Position',[10 f1.Position(4)-40 800 15],'String',[beamPhasesAngles,';  ',beamDirections], ...
        'HorizontalAlignment', 'Left','Background','w'); 


    % Figure as a 2d color map
    f2 = figure;
%     set(f2,'Position',[850,950-figHeight,800,figHeight])
    set(f2,'color','w')
    s = surf(X,Y,(avgI));
    xlabel('X Position, [$\lambda$]','interpreter','latex')
    ylabel('Y Position, [$\lambda$]','interpreter','latex')
    view(2)
    s.EdgeColor = 'none';
    colorbar
    colormap(parula)

    % Adding parameters to the top of the figure
    uicontrol(f2,'style','text','Units','Pixels',...
        'Position',[10 f2.Position(4)-20 800 15],'String',[EAmps,';  ',beamPolarizationAngles,';  '], ...
        'HorizontalAlignment', 'Left','Background','w'); 
    uicontrol(f2,'style','text','Units','Pixels',...
        'Position',[10 f2.Position(4)-40 800 15],'String',[beamPhasesAngles,';  ',beamDirections], ...
        'HorizontalAlignment', 'Left','Background','w'); 

    %%%%  This was for some 1D plots
    % f9 = figure(9);
    % set(f9,'Position',[100,900-figHeight,800,figHeight])
    % plot(X(1,:),avgI(1,:))
    % hold on
    % plot(X(1,:),avgI(11,:))
    % plot(X(1,:),avgI(21,:))
    % plot(X(1,:),avgI(31,:))
    % plot(X(1,:),avgI(41,:))
    
    
    %% Creating figures for graphs of cosine terms in the potential
    kXs = kMag*cos(th);
    kYs = kMag*sin(th);
    kVects = [transpose(kXs),transpose(kYs)];
    
    numComponentPotentials = (numWaves*(numWaves-1))/2;
    [dim1,dim2]=goodSubFigDims(numComponentPotentials);
    
    % Making larger figure
%     figHeight2 = 800;
%     figWidth = figHeight2*dim1/dim2 + 100*dim1;  % Extra 100*dim1 to account for color bar width
    
    f3 = figure();
%     set(f3,'Position',[40,950-figHeight2,figWidth,figHeight2])
    set(f3,'color','w')
    
    % Initializing parameters for the cosine terms
    subPlotLabels = cell(dim2,dim1);
    waveAmplitudes = zeros(1,numComponentPotentials);
    deltaKs = zeros(numComponentPotentials,2);
    deltaKsUnitless = zeros(numComponentPotentials,2);
    deltaPhis = zeros(1,numComponentPotentials);
    subplots = cell(1,numComponentPotentials);
    
    jj = 1; % An index for which component cosine in the potential
    
    totalFromCompSum = zeros(size(X));
    
    for step = 1:floor(numWaves/2)
        pp = 1;
        while (pp<=numWaves)&&((~(step==numWaves/2))||(pp<=numWaves/2))
            % Conditional expression:  Usually you want pp to run from 1 to
            % numWaves (this is the "(pp<=numWaves)"),  However, in the
            % case that step == numWaves/2, you will double count pairs
            % unless you stop at (pp<=numWaves/2)  (this is the "(~(step==numWaves/2))||(pp<=numWaves/2)")
            
            qq = mod(pp+step-1,numWaves)+1;  % Index of the second wave for each pair
            
            
            % Calculate parameters of the cosine term of the potential
            waveAmplitudes(jj) = EzAmp(pp)*EzAmp(qq) + Epara(pp)*Epara(qq)*cos(th(pp)-th(qq));  
            
            deltaKs(jj,:) = kVects(pp,:) - kVects(qq,:);
            deltaKsUnitless(jj,:) = deltaKs(jj,:)./(kMag);
            deltaPhis(jj) = ph(pp) - ph(qq);
            subPlotLabels{jj} = {[num2str(pp) '-' num2str(qq) ';  Depth=' num2str(2*waveAmplitudes(jj))...
                ';  $\delta\mathbf{K}$=[' num2str(round(deltaKs(jj,1)./(kMag),2)) ', ' num2str(round(deltaKs(jj,2)./(kMag),2)) '], Units of [$k_l$] ' ]
                [ '  Direction Angle:  ' num2str(round(180/pi*atan2(deltaKs(jj,2),deltaKs(jj,1))))...
                '$^{\circ}$'  '; Phase=' num2str(round((180/pi)*deltaPhis(jj),2)) '$^{\circ}$']};

            dKx = deltaKs(jj,1);
            dKy = deltaKs(jj,2);
            compWaveMatrix = waveAmplitudes(jj)*cos(dKx*X + dKy*Y - deltaPhis(jj));
            totalFromCompSum = totalFromCompSum + compWaveMatrix;

            
            % Create subplot of cosine term in potential
            subplots{jj} = subplot(dim2,dim1,jj);
%             con = contourf(X,Y,(totalFromCompSum),N);
            contourf(X,Y,(compWaveMatrix+abs(waveAmplitudes(jj))),N);
            xlabel('X Position, [$\lambda$]','interpreter','latex')
            ylabel('Y Position, [$\lambda$]','interpreter','latex')
            view(2)
            title(subPlotLabels{jj},'interpreter','latex')
            colorbar
            
            
            hold on
            
            iPts = GuideLineBoundaryIntersectionPoints(deltaKs(jj,:),deltaPhis(jj),xMin,xMax,yMin,yMax,waveAmplitudes(jj)<0);
            for rr = 1:length(iPts)
                plot([iPts(rr).x1,iPts(rr).x2],[iPts(rr).y1,iPts(rr).y2],'k','LineWidth',2)
            end
            
            % Increasing iteration index
            jj = jj+1;
            pp = pp+1;
            
%             disp('While Logical')
%             disp((pp<=numWaves)&&((~(step==numWaves/2))||(pp<=numWaves/2)))
        end
    end
    totalFromCompSum = totalFromCompSum + 0.5*sum(A.*A);

    ASorted = sort(A,'descend');
    
    climMax = 2*ASorted(1)*ASorted(2);

    for jj = 1:numComponentPotentials
        subplots{jj}.CLim = [0,climMax];
    end
    
    % Making a plot of the whole potential from the sum of component waves
    f4 = figure;
%     set(f4,'Position',[200,950-figHeight-200,800,figHeight])
    set(f4,'color','w')
    contourf(X,Y,(totalFromCompSum),N);
    xlabel('X Position, [$\lambda$]','interpreter','latex')
    ylabel('Y Position, [$\lambda$]','interpreter','latex')
    view(2)
%     s.EdgeColor = 'none';
    colorbar
    colormap(parula)
    hold on
    
    %% Drawing lines on the total plot
%     plot([0.5,1],[0.5,1],'k','LineWidth',2)
    
    
    for jj=1:numComponentPotentials
        iPts = GuideLineBoundaryIntersectionPoints(deltaKs(jj,:),deltaPhis(jj),xMin,xMax,yMin,yMax,(waveAmplitudes(jj)<0));
        for rr = 1:length(iPts)
            if ismember(jj,[1,2])
                plot([iPts(rr).x1,iPts(rr).x2],[iPts(rr).y1,iPts(rr).y2],'k','LineWidth',2)
            end
            if ismember(jj,[3,4])
                plot([iPts(rr).x1,iPts(rr).x2],[iPts(rr).y1,iPts(rr).y2],'g-','LineWidth',2)
            end
            if ismember(jj,[5,6])
                plot([iPts(rr).x1,iPts(rr).x2],[iPts(rr).y1,iPts(rr).y2],'r-','LineWidth',2)
            end
        end
        
    end
    
    

    % Adding parameters to the top of the figure
    uicontrol(f2,'style','text','Units','Pixels',...
        'Position',[10 f2.Position(4)-20 800 15],'String',[EAmps,';  ',beamPolarizationAngles,';  '], ...
        'HorizontalAlignment', 'Left','Background','w'); 
    uicontrol(f2,'style','text','Units','Pixels',...
        'Position',[10 f2.Position(4)-40 800 15],'String',[beamPhasesAngles,';  ',beamDirections], ...
        'HorizontalAlignment', 'Left','Background','w'); 
    
% %     Check Method Agreement
%     f5 = figure;
%     set(f5,'Position',[250,950-figHeight-250,800,figHeight])
%     set(f5,'color','w')
%     s = surf(X,Y,totalFromCompSum-avgI);
%     view(2)
%     s.EdgeColor = 'none';
%     colorbar
%     colormap(parula)
    

end

%% Local Functions

function [dim1,dim2] = goodSubFigDims(x)
    
    flag = 1;
    while flag

        xPrimeFacts = factor(x);

        if length(xPrimeFacts)==1
            dim1 = xPrimeFacts;
            dim2 = 1;
        else
            dim1 = xPrimeFacts(length(xPrimeFacts));
            dim2 = xPrimeFacts(length(xPrimeFacts)-1);
            for kk = fliplr(1:(length(xPrimeFacts)-2))
                if dim1 < dim2
                    dim1 = dim1*xPrimeFacts(kk);
                else
                    dim2 = dim2*xPrimeFacts(kk);
                end
            end
            if dim1 < dim2
                % Re-define dim1 and dim2 so that dim1 > dim2
                temp = dim2;
                dim2 = dim1;
                dim1 = temp;
            end
        end
        
        if dim1/dim2<2
            flag = 0;
        else
            x = x+1;
        end
    end
end


function intersectionPoints = GuideLineBoundaryIntersectionPoints(dKvec,dPh,xMin,xMax,yMin,yMax,ampIsNegative)
    % This function solves for where the cosine terms in the potential are
    % maximized.  I.e. it solves for lines along which 
    % dKvec.*[x,y] + dPh = n*2*pi OR dKvec.*[x,y] + dPh = n*2*pi + pi if
    % amp is negative
    % The output, intersectionPoints, is a struct array of the the points 
    % on the boundary of the plotted region where these lines intersect.  
    
    dPhEff = dPh + ampIsNegative*pi;   % For the case that the amplitude is negative, we run the same code but use  dPhEff = (dPh + pi) instead of dPhEff = dPh  (it is the effective dPh) 
    
    if sum(abs(dKvec))==0
        intersectionPoints = struct('x1',{},'y1',{},'x2',{},'y2',{});
        
    else
        
        boundaryCorners = [xMin,yMin;...
                        xMax,yMin;...
                        xMin,yMax;...
                        xMax,yMax];

        if abs(dKvec(2))>abs(dKvec(1))
            % Want to avoid dKvec(2)==0, and this is a way. 
            r0 = [0,dPhEff/dKvec(2)]; % r0 is just a point on the line where dKvec.*[x,y] + dPh = 0
        else
            r0 = [dPhEff/dKvec(1),0];
        end

        dispVects = zeros(4,2); % Vectors from r0 to the boundary corners
        for ii = 1:4
            dispVects(ii,:) = boundaryCorners(ii,:) - r0;
        end

        dKHat = dKvec/norm(dKvec);
        dL = 2*pi/norm(dKvec);
        
        distances = zeros(1,4);
        for ii = 1:4
            distances(ii) =  sum(dispVects(ii,:).*dKHat);  % Distances from the corners to the line corresponding to dKvec.*[x,y] + dPh = 0
        end
        nMin = ceil( min(distances) / dL );  % Lowest value of n that gives a line dKvec.*[x,y] + dPh = n*2*pi that crosses through the boundary
        nMax = floor( max(distances) / dL ); % Highest value of n that gives a line dKvec.*[x,y] + dPh = n*2*pi that crosses through the boundary

        nLines = nMax-nMin+1;
        intersectionPoints = struct('dKvec',cell(1,nLines),'dPh',cell(1,nLines),'x1',cell(1,nLines),'y1',cell(1,nLines),'x2',cell(1,nLines),'y2',cell(1,nLines));
    
        [intersectionPoints.dKvec] = deal(dKvec);
        [intersectionPoints.dPh] = deal(dPh);
        
        for n = nMin:nMax
            index = 1;
            % Check intersection with x=xMin
            xMinIntersect = (-(dKvec(1)/dKvec(2))*xMin + (dPhEff + n*2*pi)/dKvec(2));
            if ((yMin <= xMinIntersect)&&(xMinIntersect< yMax))
                intersectionPoints(n-nMin+1).x1 = xMin;
                intersectionPoints(n-nMin+1).y1 = xMinIntersect;
                index = index + 1;
            end
            
            % Check intersection with x=xMax
            xMaxIntersect = (-(dKvec(1)/dKvec(2))*xMax + (dPhEff + n*2*pi)/dKvec(2));
            if ((yMin <= xMaxIntersect)&&(xMaxIntersect < yMax))
                if index==1
                    intersectionPoints(n-nMin+1).x1 = xMax;
                    intersectionPoints(n-nMin+1).y1 = xMaxIntersect;
                elseif index==2
                    intersectionPoints(n-nMin+1).x2 = xMax;
                    intersectionPoints(n-nMin+1).y2 = xMaxIntersect;
                end
                index = index + 1;
            end
            
            % Check intersection with y=yMin
            yMinIntersect = (-(dKvec(2)/dKvec(1))*yMin + (dPhEff + n*2*pi)/dKvec(1));
            if ((xMin < yMinIntersect)&&(yMinIntersect <= xMax))
                if index==1
                    intersectionPoints(n-nMin+1).x1 = yMinIntersect;
                    intersectionPoints(n-nMin+1).y1 = yMin;
                elseif index==2
                    intersectionPoints(n-nMin+1).x2 = yMinIntersect;
                    intersectionPoints(n-nMin+1).y2 = yMin;
                end
                
                index = index + 1;
            end
            
            % Check intersection with y=yMax
            yMaxIntersect = (-(dKvec(2)/dKvec(1))*yMax + (dPhEff + n*2*pi)/dKvec(1));
            if ((xMin < yMaxIntersect)&&(yMaxIntersect <= xMax))
                if index==1
                    intersectionPoints(n-nMin+1).x1 = yMaxIntersect;
                    intersectionPoints(n-nMin+1).y1 = yMax;
                elseif index==2
                    intersectionPoints(n-nMin+1).x2 = yMaxIntersect;
                    intersectionPoints(n-nMin+1).y2 = yMax;
                end
            end      
        end
    end

end
