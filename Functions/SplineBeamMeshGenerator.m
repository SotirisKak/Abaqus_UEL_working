function [BeamMesh]= SplineBeamMeshGenerator(BeamAxes, BeamOrder, hbeam)
% Generates Beam Mesh depending on anything dependable
n_beams = size(BeamAxes.cs,2);

%Export Gauss points for length calculations (quadratic elems)
deg_threshold = 40; % Iterative straighten high distorted elements
nGP = 12;
[xis, wGP] = lgwt(nGP,-1,1);
[Phi, dPhi] = ShapeFnc(xis,  BeamOrder);
rdp = BeamAxes.rdp;
% Loop over beams
iseg = 0;
for i = 1:n_beams
    cs = BeamAxes.cs(i);
    seg_mat = BeamAxes.segment_mtrx{i};
    n_seg = size(seg_mat,1);
    
    for j = 1:n_seg
        iseg = iseg+1;
        lb = seg_mat(j,4); % Beam length
        nh = max([round(lb/hbeam),1]); % Beam elements
        le = lb/nh; % Length element
        IDcb = seg_mat(j,3);% RVE

        xa = seg_mat(j,1); xb = seg_mat(j,2);
        Dx = (xb-xa)/nh/100;
        x_dis = linspace(xa,xb,round((xb-xa)/Dx));

        % Initialize
        pxyz = zeros(length(x_dis),3);
        s_loc = zeros(length(x_dis),1)';
        pxyz(1,:) = fnval(cs,xa)'-rdp(IDcb,:);
        s_loc(1) = 0;
        for ix = 2:length(x_dis)
            pxyz(ix,:) = fnval(cs,x_dis(ix))'-rdp(IDcb,:);
            s_loc(1,ix) = arclength(pxyz(1:ix,1),pxyz(1:ix,2),pxyz(1:ix,3));
        end

      
        % Loop over elements
        SParam = zeros(nh,BeamOrder+1);
        ArcLength = 0;
        FE_ArcLength = 0;
        FESParam = zeros(nh,BeamOrder+1);
        x = interp1(s_loc,x_dis,ArcLength);
        Nodes = fnval(cs,x)'-rdp(IDcb,:)-[0,0,1];
        Nodes(Nodes<-1) = -1;
        Nodes(Nodes>1) = 1;
        
        Connectivity = [];
        
        for ii = 1:nh
        
            if BeamOrder == 1 || BeamOrder==3%Linear Elements OR EulerBernoulli
                            
                SParam(ii,:) = [ArcLength,ArcLength + le];
                ArcLength = ArcLength + le;
                x = interp1(s_loc,x_dis,ArcLength);
                rc = (fnval(cs,x)'-rdp(IDcb,:)-[0,0,1])';
                rc(rc<-1) = -1;
                rc(rc>1) = 1;
                
                % Update discretized arc length
                FEle = norm(Nodes(end,:)-rc');
                FESParam(ii,:) = [FE_ArcLength,FE_ArcLength + FEle];
                FE_ArcLength = FE_ArcLength + FEle;
                
                % Update nodes list and connectivity        
                Nodes = [Nodes; rc'];
                Connectivity = [Connectivity; ii, ii+1];
                
            elseif BeamOrder == 2 %Quadratic
                SParam(ii,:) = [ArcLength,ArcLength+le/2,ArcLength + le];
                MidArcLength =  ArcLength + le/2;
                ArcLength = ArcLength + le;
      
                x = interp1(s_loc,x_dis,MidArcLength);
                rm = (fnval(cs,x)'-rdp(IDcb,:)-[0,0,1])';
                rm(rm<-1) = -1;
                rm(rm>1) = 1;
                
                if ArcLength<=s_loc(end)
                    x = interp1(s_loc,x_dis,ArcLength);
                else
                    x = xb;
                end

                rc = (fnval(cs,x)'-rdp(IDcb,:)-[0,0,1])';
                rc(rc<-1) = -1;
                rc(rc>1) = 1;
                
                % Straighten up distorted elements
                ra = Nodes(end,:)';
                rm_temp = rm;
                rb = rc;
                
                u = rm_temp-ra; v= rb-rm_temp;
                an_theta = rad2deg(atan2(norm(cross(u,v)),dot(u,v)));
                
                while an_theta>deg_threshold
                    
                    rm_temp = (rm_temp+(ra+rb)/2)/2;
                    u = rm_temp-ra; v= rb-rm_temp;
                    an_theta = rad2deg(atan2(norm(cross(u,v)),dot(u,v)));
%                     figure(1)
%                     hold on
%                     plot3([ra(1),rm_temp(1),rb(1)],[ra(2),rm_temp(2),rb(2)],[ra(3),rm_temp(3),rb(3)])
                    
                end
                rm = rm_temp; % Updated, non-distorted element
                
                
                % Update discretized arc length
                locNodes = [Nodes(end,:); rm'; rc']';      
                fd = (dPhi*locNodes');
                JA2 = (sum(fd.^2,2)).^(1/2);
                FEle = sum(JA2.*wGP);
                
                FESParam(ii,:) = [FE_ArcLength,FE_ArcLength+FEle/2,FE_ArcLength + FEle];
                FE_ArcLength = FE_ArcLength + FEle;
                
                
                % Keep Nodes and Connectivity
                Nodes = [Nodes; rm'; rc'];
                istart = (ii-1)*2+1;
                Connectivity = [Connectivity; istart, istart+1, istart+2];

                if sum(sum(isnan(Nodes)))~=0
                    stop
                end
            end
            
        end
        % HERE FOR BEAM NETWORKS. Will need to figure out sth for common nodes for
        % the SParam vector, since if I merge the nodes then will also have to
        % include indices as well there.
        BeamMesh(iseg).Nodes = Nodes;
        BeamMesh(iseg).SParam = FESParam;
        BeamMesh(iseg).PhysicalSParam = SParam;
        BeamMesh(iseg).Connectivity = Connectivity;
        BeamMesh(iseg).nB = size(Nodes,1);
        BeamMesh(iseg).Order = BeamOrder;
    end

end





    
    
    
    
    
    
    
    
