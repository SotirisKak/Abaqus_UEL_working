function L2 = ufemError2(gp,w,RefU,RefNodes,RefCon,CurU,CurNodes,CurCon, GlobalNodes, GlobalCon)

    % Calculate centroids of bricks
    for i = 1:size(RefCon,1)
        temp = RefNodes(RefCon(i,:),:);
        refSCentroids(i,:) = mean(temp);
    end
    refnSel = size(RefCon,1);
    
    for i = 1:size(CurCon,1)
        temp = CurNodes(CurCon(i,:),:);
        SCentroids(i,:) = mean(temp);
    end
    nSel = size(CurCon,1);

    % Initialize
    L2 = 0;
   

     % Progress bar
%         waitbar(i/size(RefCon,1),fprogressbar)

    % Export Integration Domain
    rk = GlobalNodes(GlobalCon,:); % Global c
   
    for j = 1:length(w)
        zetas = gp(j,:);
        [Nks, dNks] = Hex8ShapeFnc(zetas(1),zetas(2),zetas(3));

        %Coordinates of GP in the physical domain
        xref = Nks'*rk;
        
        %% REFERENCE DISPLACEMENT
        % Find 6 closest centroids
        if refnSel>6
            [IDX, ~] = knnsearch(refSCentroids,xref,'K',6);
        else
            IDX = 1:refnSel;
        end

        ID0 = PointInBrick(xref', RefCon(IDX,:), RefNodes);
        if isempty(ID0)
           error('SKATA')
        end

        InBrick =IDX(ID0);
        BrickIDs = RefCon(InBrick(1),:);

        rk_ref = RefNodes(BrickIDs,:); % Global coordinates 
        [xi_ref,zeta_ref,eta_ref] = InverseMapping(xref', rk_ref, 1);
        Nks_ref = Hex8ShapeFnc(xi_ref, zeta_ref, eta_ref);

        xtestref = Nks_ref'*rk_ref;
        uref = Nks_ref'*RefU(BrickIDs,2:end); % displacement

        %% NOW DEAL WITH THE COMPARISON SOLUTION
        % Find 6 closest centroids
        if nSel>6
            [IDX, ~] = knnsearch(SCentroids,xref,'K',6);
        else
            IDX = 1:nSel;
        end

        ID0 = PointInBrick(xref', CurCon(IDX,:), CurNodes);
        if isempty(ID0)
           error('SKATA')
        end

        InBrick =IDX(ID0);
        BrickIDs = CurCon(InBrick(1),:);

        rk_cur = CurNodes(BrickIDs,:); % Global coordinates 
        [xi_cur,zeta_cur,eta_cur] = InverseMapping(xref', rk_cur, 1);
        Nks_cur = Hex8ShapeFnc(xi_cur, zeta_cur, eta_cur);

        xtestcur = Nks_cur'*rk_cur;
        ucomp = Nks_cur'*CurU(BrickIDs,2:end); % displacement

        % Perform the integration
        Jdet = det(rk'*dNks);
        L2 = L2 + w(j)*Jdet* norm(ucomp-uref)^2;
%          L2 = L2 + w(j)*Jdet;


    end
       
            
        