function [] = D01_WriteAbaqusInput_RPL_PBC(WorkDir, JobNum, SolidMesh, BeamMesh, PenaltyConst, BeamRadius, param_beam, X_disp, nSteps, Matl, param_solid, beam_nu, CoupleElemCon, nmax, nT, SolverOptions,RPL)


filename = [WorkDir, '/Job',num2str(JobNum,'%.4d'),'/Model',num2str(JobNum,'%.4d'),RPL,'.inp'];
fid = fopen(filename,'w');

%% GEOMETRY SECTION
BeamOrder = BeamMesh.Order;
nS = size(SolidMesh.Nodes,1);
nB = size(BeamMesh.Nodes,1);
nh = size(BeamMesh.Connectivity,1);

Nodes = [SolidMesh.Nodes; BeamMesh.Nodes];
ntotal = nS+nB;
GeometryData = [(1:ntotal)',Nodes];

% find free end bricks
[YCs,S1FaceBricks]=sort(SolidMesh.Centroids(:,3));
S1FaceBricks = S1FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S2FaceBricks]=sort(SolidMesh.Centroids(:,3),'descend');
S2FaceBricks = S2FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S6FaceBricks]=sort(SolidMesh.Centroids(:,1));
S6FaceBricks = S6FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S4FaceBricks]=sort(SolidMesh.Centroids(:,1),'descend');
S4FaceBricks = S4FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S3FaceBricks]=sort(SolidMesh.Centroids(:,2));
S3FaceBricks = S3FaceBricks(1:length(YCs)/length(unique(YCs)));

[YCs,S5FaceBricks]=sort(SolidMesh.Centroids(:,2),'descend');
S5FaceBricks = S5FaceBricks(1:length(YCs)/length(unique(YCs)));



W=max(Nodes(:,1))-min(Nodes(:,1));
D=max(Nodes(:,2))-min(Nodes(:,2));
H=max(Nodes(:,3))-min(Nodes(:,3));

fprintf(fid, '*NODE \n');
fprintf(fid, '\t %i, %2.14e, %2.14e, %2.14e \n', GeometryData');

fprintf(fid, '*NODE, NSET=DEF\n');
fprintf(fid, '\t%d,0,0,0\n',ntotal+1);
fprintf(fid, '\t%d,0,0,0\n',ntotal+2);
fprintf(fid, '\t%d,0,0,0\n',ntotal+3);
fprintf(fid, '\t%d,0,0,0\n',ntotal+4);
fprintf(fid, '\t%d,0,0,0\n',ntotal+5);
fprintf(fid, '\t%d,0,0,0\n',ntotal+6);

%% PBC NODE SETS
% Faces
PX = SolidMesh.posx; NX = SolidMesh.negx;
PY = SolidMesh.posy; NY = SolidMesh.negy;
PZ = SolidMesh.posz; NZ = SolidMesh.negz;

% Corners
C1 = intersect(intersect(PX,PY),PZ); C5 = intersect(intersect(PX,NY),PZ);
C2 = intersect(intersect(PX,PY),NZ); C6 = intersect(intersect(PX,NY),NZ);
C3 = intersect(intersect(NX,PY),NZ); C7 = intersect(intersect(NX,NY),NZ);
C4 = intersect(intersect(NX,PY),PZ); C8 = intersect(intersect(NX,NY),PZ);

% Edges
E1 = setxor(intersect(PX,PY),[C1,C2]); E5 = setxor(intersect(PX,NY),[C5,C6]);
E2 = setxor(intersect(NZ,PY),[C2,C3]); E6 = setxor(intersect(NZ,NY),[C6,C7]);
E3 = setxor(intersect(NX,PY),[C3,C4]); E7 = setxor(intersect(NX,NY),[C7,C8]);
E4 = setxor(intersect(PZ,PY),[C4,C1]); E8 = setxor(intersect(PZ,NY),[C8,C5]);

E9 = setxor(intersect(PX,PZ),[C1,C5]);
E10 = setxor(intersect(PX,NZ),[C2,C6]);
E11 = setxor(intersect(NX,NZ),[C3,C7]);
E12 = setxor(intersect(NX,PZ),[C4,C8]);

% Faces
F1 = setxor(PX,[C1,E1,C2,E10,C6,E5,C5,E9]);
F2 = setxor(NZ,[C2,E2,C3,E11,C7,E6,C6,E10]);
F3 = setxor(NX,[C3,E3,C4,E12,C8,E7,C7,E11]);
F4 = setxor(PZ,[C4,E4,C1,E9,C5,E8,C8,E12]);

F5 = setxor(PY,[C1,E1,C2,E2,C3,E3,C4,E4]);
F6 = setxor(NY,[C5,E5,C6,E6,C7,E7,C8,E8]);

% BC Sets
F1_BC = [F1,E1,E9,E5,E10,C1,C5,C6,C2];
F2_BC = [F2,E2,E10,E6,E11,C2,C6,C7,C3];
F3_BC = [F3,E3,E12,E7,E11,C3,C4,C8,C7];
F4_BC = [F4,E8,E9,E4,E12,C4,C8,C5,C1];
F5_BC = [F5,E4,E1,E2,E3,C1,C2,C3,C4];
F6_BC = [F6,E5,E6,E7,E8,C5,C6,C7,C8];

% Align Sets (register corresponding points)
F5 = nodeRegister(F5,F6,SolidMesh.Nodes);
F1 = nodeRegister(F1,F3,SolidMesh.Nodes);
F2 = nodeRegister(F2,F4,SolidMesh.Nodes);

F5_BC = nodeRegister(F5_BC,F6_BC,SolidMesh.Nodes);
F1_BC = nodeRegister(F1_BC,F3_BC,SolidMesh.Nodes);
F2_BC = nodeRegister(F2_BC,F4_BC,SolidMesh.Nodes);

E3 = nodeRegister(E3,E1,SolidMesh.Nodes);
E7 = nodeRegister(E7,E3,SolidMesh.Nodes);
E5 = nodeRegister(E5,E7,SolidMesh.Nodes);

E12 = nodeRegister(E12,E9,SolidMesh.Nodes);
E11 = nodeRegister(E11,E12,SolidMesh.Nodes);
E10 = nodeRegister(E10,E11,SolidMesh.Nodes);

E8 = nodeRegister(E8,E4,SolidMesh.Nodes);
E6 = nodeRegister(E6,E8,SolidMesh.Nodes);
E2 = nodeRegister(E2,E6,SolidMesh.Nodes);


% WRITE NODE SETS
writeNSET(fid, F5, 'TopS'); writeNSET(fid, F6, 'BotS');
writeNSET(fid, F1, 'FrontS'); writeNSET(fid, F3, 'BackS');
writeNSET(fid, F4, 'LeftS'); writeNSET(fid, F2, 'RightS')

writeNSET(fid, F5_BC, 'Top_BC'); writeNSET(fid, F6_BC, 'Bot_BC');
writeNSET(fid, F1_BC, 'Front_BC'); writeNSET(fid, F3_BC, 'Back_BC');
writeNSET(fid, F4_BC, 'Left_BC'); writeNSET(fid, F2_BC, 'Right_BC');

writeNSET(fid, E1, 'FTedge'); writeNSET(fid, E3, 'BTedge'); 
writeNSET(fid, E7, 'BBedge'); writeNSET(fid, E5, 'FBedge');

writeNSET(fid, E9, 'FLedge'); writeNSET(fid, E12, 'BLedge'); 
writeNSET(fid, E11,'BRedge'); writeNSET(fid, E10, 'FRedge');

writeNSET(fid, E4, 'LTedge'); writeNSET(fid, E8, 'LBedge'); 
writeNSET(fid, E6, 'RBedge'); writeNSET(fid, E2, 'RTedge');

writeNSET(fid, C1, 'C1'); writeNSET(fid, C4, 'C2'); 
writeNSET(fid, C3, 'C3'); writeNSET(fid, C2, 'C4');
writeNSET(fid, C5, 'C5'); writeNSET(fid, C8, 'C6'); 
writeNSET(fid, C7, 'C7'); writeNSET(fid, C6, 'C8');




%% Ref Node Sets
fprintf(fid, '*NSET, NSET=FixRefNode\n');
fprintf(fid, '%i \n', ntotal+1);
fprintf(fid, '*NSET, NSET=FreeRefNode\n');
fprintf(fid, '%i \n', ntotal+2);



n_count = 0;
%% ElEMENT SECTION
% Declare UEL
fprintf(fid, '*USER ELEMENT, TYPE=U001, NODES=%i, PROPERTIES=2, COORDINATES=3, UNSYMM \n', nmax.A);
fprintf(fid, '1,2,3 \n'); 
% Declare UEL
fprintf(fid, '*USER ELEMENT, TYPE=U002, NODES=%i, PROPERTIES=2, COORDINATES=3, UNSYMM \n', nmax.B);
fprintf(fid, '1,2,3 \n');
% Declare UEL
fprintf(fid, '*USER ELEMENT, TYPE=U003, NODES=%i, PROPERTIES=2, COORDINATES=3, UNSYMM \n', nmax.C);
fprintf(fid, '1,2,3 \n');


% Bricks
if strcmp(Matl, 'HYPERFOAM')
    fprintf(fid, '*ELEMENT, TYPE=C3D8R, ELSET=BrickElems \n');
else
    fprintf(fid, '*ELEMENT, TYPE=C3D8RH, ELSET=BrickElems \n');
end

SCon = SolidMesh.Connectivity;
nSolEl = size(SCon,1); 
SolidConnectivity = [(1:nSolEl)',SCon];
fprintf(fid, ['%i', repmat(', %i',1,8), ' \n'], SolidConnectivity');

% Beams
BCon = BeamMesh.Connectivity + nS; % add nS for global connectivity
nBeamEl = size(BCon,1); 
nTotEl = nBeamEl + nSolEl;
BeamConnectivity = [(nSolEl+1:nTotEl)',BCon];
% fprintf(fid, '*ELEMENT, TYPE=B3%d , ELSET=BeamElems \n', BeamOrder);
% if BeamOrder == 1 || BeamOrder == 3
%     fprintf(fid, ['%i', repmat(', %i',1,2), ' \n'], BeamConnectivity');
% else
%     fprintf(fid, ['%i', repmat(', %i',1,3), ' \n'], BeamConnectivity');
% end
for i = 1:nh
    fprintf(fid, '*ELEMENT, TYPE=B3%dH, ELSET=BeamElem%d \n', [BeamOrder,i]);
    if BeamOrder == 1 || BeamOrder == 3
        fprintf(fid, ['%i', repmat(', %i',1,2), ' \n'], BeamConnectivity(i,:)');
    else
        fprintf(fid, ['%i', repmat(', %i',1,3), ' \n'], BeamConnectivity(i,:)');
    end
end


                
    

% Define normals nq
fprintf(fid, '*NORMAL \n');
for i = 1:nh
    fprintf(fid, ['%d, %d, %.6f, %.6f, %.6f \n'], [BeamConnectivity(i,1),BeamConnectivity(i,2),BeamMesh.n2(i,:)]);
%     % Assign the same normal to the end point as the beginning of the next
%     % element
%     if i<nh
%         if BeamConnectivity(i,end)==BeamConnectivity(i+1,1)
%             fprintf(fid, ['%d, %d, %.6f, %.6f, %.6f \n'], [BeamConnectivity(i,1),BeamConnectivity(i,end),BeamMesh.n2(i+1,:)]);
%         end
%     end
            
end

% UELA
ncplA = nmax.A;
for i = 1:size(CoupleElemCon.A,1)
    fprintf(fid, '*ELEMENT, TYPE=U001 , ELSET=CplEl%d \n',i);    
    UEL_Con = [nTotEl+i, CoupleElemCon.A(i,:)];
    
    if (ncplA+1) < 8
        str_line = ['%i', repmat(', %i',1,nmax.A), ' \n'];
        fprintf(fid, str_line, UEL_Con);
    else
        n8rows = fix((ncplA+1)/8);
        nrest = mod((ncplA+1),8);
        str_line = ['%i', repmat(', %i',1,7), ', \n'];
        fprintf(fid, str_line, UEL_Con(1:8*n8rows));
        if nrest>0
            str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
            fprintf(fid, str_line, UEL_Con(8*n8rows+1:end));
        end
    end
end

% UELB
ncplB = nmax.B;
for i = 1:size(CoupleElemCon.B,1)
    fprintf(fid, '*ELEMENT, TYPE=U002 , ELSET=CplEl%d \n',i+size(CoupleElemCon.A,1));    
    UEL_Con = [nTotEl+i+size(CoupleElemCon.A,1), CoupleElemCon.B(i,:)];
    
    if (ncplB+1) < 8
        str_line = ['%i', repmat(', %i',1,nmax.B), ' \n'];
        fprintf(fid, str_line, UEL_Con);
    else
        n8rows = fix((ncplB+1)/8);
        nrest = mod((ncplB+1),8);
        str_line = ['%i', repmat(', %i',1,7), ', \n'];
        fprintf(fid, str_line, UEL_Con(1:8*n8rows));
        if nrest>0
            str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
            fprintf(fid, str_line, UEL_Con(8*n8rows+1:end));
        end
    end
end

% UELC
ncplC = nmax.C;
for i = 1:size(CoupleElemCon.C,1)
    fprintf(fid, '*ELEMENT, TYPE=U003 , ELSET=CplEl%d \n',i+size(CoupleElemCon.A,1)+size(CoupleElemCon.B,1));    
    UEL_Con = [nTotEl+i+size(CoupleElemCon.A,1)+size(CoupleElemCon.B,1), CoupleElemCon.C(i,:)];
    
    if (ncplC+1) < 8
        str_line = ['%i', repmat(', %i',1,nmax.C), ' \n'];
        fprintf(fid, str_line, UEL_Con);
    else
        n8rows = fix((ncplC+1)/8);
        nrest = mod((ncplC+1),8);
        str_line = ['%i', repmat(', %i',1,7), ', \n'];
        fprintf(fid, str_line, UEL_Con(1:8*n8rows));
        if nrest>0
            str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
            fprintf(fid, str_line, UEL_Con(8*n8rows+1:end));
        end
    end
end

% if crosslinks.merge
%     
%     %fprintf(fid, '*ELEMENT, TYPE=CONN3D2, ELSET=Crosslinks \n');
%     
%     
%     % Loop through merged nodes groups
%     n_count = 0;
%     for i = 1:size(BeamMesh.Nodes2Merge,1)
%         
%         n_loc_merge = BeamMesh.Nodes2Merge{i, 2};
%         
%         for iii = 1:crosslinks.dof
%             fprintf(fid, '*ELEMENT, TYPE=SPRING2, ELSET=Crosslinks%d\n',iii);
%         
%             for ii = 1:n_loc_merge
% 
%                 n_count = n_count + 1;
%                 n_loc_nodes = BeamMesh.Nodes2Merge{i, 1} + nS;
% 
%                 if ii<n_loc_merge
% 
%                     fprintf(fid, '%d, %d, %d\n',[nTotEl+nh+n_count,n_loc_nodes(ii),n_loc_nodes(ii+1)]);
% 
%                 elseif ii>=n_loc_merge && n_loc_merge>2
% 
%                     fprintf(fid, '%d, %d, %d\n',[nTotEl+nh+n_count,n_loc_nodes(ii),n_loc_nodes(1)]);
% 
%                 end
%             end
%         end
%     end
% end

for i = 1:size(CoupleElemCon.A,1)
    fprintf(fid, '*UEL PROPERTY, ELSET=CplEl%d \n',i);
    fprintf(fid, '%2.3E, %d \n', [PenaltyConst, i]);
end

for i = 1:size(CoupleElemCon.B,1)
    fprintf(fid, '*UEL PROPERTY, ELSET=CplEl%d \n',i+size(CoupleElemCon.A,1));
    fprintf(fid, '%2.3E, %d \n', [PenaltyConst, i]);
end

for i = 1:size(CoupleElemCon.C,1)
    fprintf(fid, '*UEL PROPERTY, ELSET=CplEl%d \n',i+size(CoupleElemCon.A,1)+size(CoupleElemCon.B,1));
    fprintf(fid, '%2.3E, %d \n', [PenaltyConst, i]);
end


%% SECTIONS
fprintf(fid, '** \n');
fprintf(fid, '** SECTIONS \n');
fprintf(fid, '** \n');
if strcmp(Matl, 'NHK') || strcmp(Matl, 'OGD') || strcmp(Matl, 'HYPERFOAM')
    fprintf(fid, '*SOLID SECTION, ELSET=BrickElems, material=SolidMatl \n');
    fprintf(fid, ',\n');
elseif strcmp(Matl, 'SVK')
    fprintf(fid, '*SOLID SECTION, ELSET=BrickElems, material=UANISO_SVK, orientation=DefaultOri \n');
    fprintf(fid, ',\n');
    fprintf(fid, '*Hourglass stiffness\n');
    fprintf(fid, '1.0,1.0,1.0\n');
else
    error('ERROR... MATERIAL NOT SUPPORTED (YET)')
end

for i = 1:nh
    fprintf(fid, '*BEAM SECTION, ELSET=BeamElem%d, material=BeamMatl, section=CIRC \n',i);
    fprintf(fid, '%2.7e \n', BeamRadius);
    % fprintf(fid, '%2.14e, %2.14e, %2.14e \n', Na);
    fprintf(fid, '%.6f, %.6f, %.6f \n',BeamMesh.n1(i,:));
end

% % Connector Section
% if crosslinks.merge
%     
% %     fprintf(fid, '*CONNECTOR SECTION, ELSET=Crosslinks \n');
% %     fprintf(fid, [crosslinks.type, '\n']);
% 
%     for iii = 1:crosslinks.dof
%     
%         fprintf(fid, '*SPRING, ELSET=Crosslinks%d \n',iii);
%         fprintf(fid, '%d, %d\n',[iii,iii]);
%         fprintf(fid, '%2.6E \n', crosslinks.stiff(iii));
%     
%     end
%     
% end

%% Constraints

fprintf(fid, '** \n');
fprintf(fid, '** Constraints \n');
fprintf(fid, '** \n');
if strcmp(RPL, 'UEx') || strcmp(RPL, 'UEy') || strcmp(RPL, 'UEz')
    writeCONSTRAINT_UE(fid,'TopS', 'BotS', 5, 1, -1, [0,-1,0],ntotal)
    writeCONSTRAINT_UE(fid,'FrontS', 'BackS', 4, 1, -1, [-1,0,0],ntotal)
    writeCONSTRAINT_UE(fid,'LeftS', 'RightS', 6, 1, -1, [0,0,-1],ntotal)
    
    writeCONSTRAINT_UE(fid,'FTedge', 'BTedge', 4, 1, -1, [-1,0,0],ntotal)
    writeCONSTRAINT_UE(fid,'BTedge', 'BBedge', 5, 1, -1, [0,-1,0],ntotal)
    writeCONSTRAINT_UE(fid,'BBedge', 'FBedge', 4, 1, -1, [1,0,0],ntotal)
    
    writeCONSTRAINT_UE(fid,'FLedge', 'BLedge', 4, 1, -1, [-1,0,0],ntotal)
    writeCONSTRAINT_UE(fid,'BLedge', 'BRedge', 6, 1, -1, [0,0,-1],ntotal)
    writeCONSTRAINT_UE(fid,'BRedge', 'FRedge', 4, 1, -1, [1,0,0],ntotal)
    
    writeCONSTRAINT_UE(fid,'LTedge', 'LBedge', 5, 1, -1, [0,-1,0],ntotal)
    writeCONSTRAINT_UE(fid,'LBedge', 'RBedge', 6, 1, -1, [0,0,-1],ntotal)
    writeCONSTRAINT_UE(fid,'RBedge', 'RTedge', 5, 1, -1, [0,1,0],ntotal)
    
    writeCONSTRAINT_UE(fid,'C6', 'C2', 5, 1, -1, [0,1,0],ntotal)
    writeCONSTRAINT_UE(fid,'C2', 'C3', 6, 1, -1, [0,0,-1],ntotal)
    writeCONSTRAINT_UE(fid,'C3', 'C4', 4, 1, -1, [1,0,0],ntotal)
    writeCONSTRAINT_UE(fid,'C4', 'C8', 5, 1, -1, [0,-1,0],ntotal)
    writeCONSTRAINT_UE(fid,'C8', 'C5', 6, 1, -1, [0,0,1],ntotal)
    writeCONSTRAINT_UE(fid,'C5', 'C1', 5, 1, -1, [0,1,0],ntotal)
    writeCONSTRAINT_UE(fid,'C1', 'C7', [4,5,6], 1, -1, [-1,-1,-1],ntotal)
elseif strcmp(RPL, 'SSxy') || strcmp(RPL, 'SSyz') || strcmp(RPL, 'SSxz')
    
    nBC_tot = length([F5_BC,F6_BC,F1_BC,F2_BC,F3_BC,F4_BC,E1,E2,E3,E4,E5,E6,E7,E9,E10,E11,E12,C1,C2,C3,C4,C5,C6,C7,C8])*3;
    fprintf(fid, '*Equation\n');
    fprintf(fid, '%d\n',9087);
    writeCONSTRAINT_UE(fid,'Top_BC', 'Bot_BC', [4,1,6], 1, -1, [-1,-1,-1],ntotal,F5_BC,F6_BC)
    writeCONSTRAINT_UE(fid,'Left_BC', 'Right_BC', [5,6,2], 1, -1, [-1,-1,-1],ntotal,F4_BC,F2_BC)
    writeCONSTRAINT_UE(fid,'Front_BC', 'Back_BC', [3,4,5], 1, -1, [-1,-1,-1],ntotal,F1_BC,F3_BC)
    
    writeCONSTRAINT_UE(fid,'FTedge', 'BTedge', [3,4,5], 1, -1, [-1,-1,-1],ntotal,E1,E3)
    writeCONSTRAINT_UE(fid,'BTedge', 'BBedge', [4,1,6], 1, -1, [-1,-1,-1],ntotal,E3,E7)
    writeCONSTRAINT_UE(fid,'BBedge', 'FBedge', [3,4,5], 1, -1, [1,1,1],ntotal,E7,E5)
    
    writeCONSTRAINT_UE(fid,'FLedge', 'BLedge', [3,4,5], 1, -1, [-1,-1,-1],ntotal,E9,E12)
    writeCONSTRAINT_UE(fid,'BLedge', 'BRedge', [5,6,2], 1, -1, [-1,-1,-1],ntotal,E12,E11)
    writeCONSTRAINT_UE(fid,'BRedge', 'FRedge', [3,4,5], 1, -1, [1,1,1],ntotal,E11,E10)
    
    writeCONSTRAINT_UE(fid,'LTedge', 'LBedge', [3,1,6], 1, -1, [-1,-1,-1],ntotal,E4,E8)
    writeCONSTRAINT_UE(fid,'LBedge', 'RBedge', [5,6,2], 1, -1, [-1,-1,-1],ntotal,E8,E6)
    writeCONSTRAINT_UE(fid,'RBedge', 'RTedge', [4,1,6], 1, -1, [1,1,1],ntotal,E6,E2)
    
    writeCONSTRAINT_UE(fid,'C6', 'C2', [4,1,6], 1, -1, [1,1,1],ntotal,C8,C4)
    writeCONSTRAINT_UE(fid,'C2', 'C3', [5,6,2], 1, -1, [-1,-1,-1],ntotal,C4,C3)
    writeCONSTRAINT_UE(fid,'C3', 'C4', [3,4,5], 1, -1, [1,1,1],ntotal,C3,C2)
    writeCONSTRAINT_UE(fid,'C4', 'C8', [4,1,6], 1, -1, [-1,-1,-1],ntotal,C2,C6)
    writeCONSTRAINT_UE(fid,'C8', 'C5', [5,6,2], 1, -1, [1,1,1],ntotal,C6,C5)
    writeCONSTRAINT_UE(fid,'C5', 'C1', [4,1,6], 1, -1, [1,1,1],ntotal,C5,C1)
    
    writeCONSTRAINT_SS(fid,'C1', 'C7', [3,1,2],[4,4,5],[5,6,6], -1,-1,-1,ntotal,C1,C7)
    
    
end

%% Materials
fprintf(fid, '** \n');
fprintf(fid, '** Materials \n');
fprintf(fid, '** \n');

if strcmp(Matl, 'NHK')%Solid Neo Hookean
    
    fprintf(fid, '*Material, name=SolidMatl \n');
    fprintf(fid, '*Hyperelastic, neo hooke \n');
    fprintf(fid, '%2.3e, %2.3e \n', param_solid);
    
elseif strcmp(Matl, 'HYPERFOAM')%Solid Neo Hookean
    
    fprintf(fid, '*Material, name=SolidMatl \n');
    fprintf(fid, '*HYPERFOAM, N=1 \n');
    fprintf(fid, '%2.3e, 2.0, %2.3e \n', param_solid);
    
elseif strcmp(Matl, 'OGD')% Ogden Incompressible
    
    fprintf(fid, '*Material, name=SolidMatl \n');
    fprintf(fid, '*Hyperelastic, OGDEN, N=1 \n');
    fprintf(fid, '%2.3e, %2.3e, %2.3e \n', param_solid);
    
elseif strcmp(Matl, 'MNRIV')% Ogden Incompressible

    fprintf(fid, '*Material, name=SolidMatl \n');
    fprintf(fid, '*Hyperelastic, MOONEY-RIVLIN \n');
    fprintf(fid, '%2.3e, %2.3e, %2.3e \n', param_solid);
    
elseif strcmp(Matl, 'SVK') %Solid SVK
    
    fprintf(fid, '*orientation, name=DefaultOri, local directions=3 \n');
    fprintf(fid, '1.0,0.0,0.0,0.0,1.0,0.0 \n');
    fprintf(fid, '3,0.0 \n');
    fprintf(fid, '1,0,0 \n');
    fprintf(fid, '0,1,0 \n');
    fprintf(fid, '0,0,1 \n');
    
    fprintf(fid, '*Material, name=UANISO_SVK \n');
    fprintf(fid, '*density \n');
    fprintf(fid, '1e-12 \n');
    
    fprintf(fid, '*ANISOTROPIC HYPERELASTIC, USER, FORMULATION=STRAIN, TYPE=COMPRESSIBLE, PROPERTIES=9 \n');
%     fprintf(fid, '%2.14e, %2.14e, %2.14e, %2.14e, %2.14e, %2.14e, %2.14e, %2.14e, %2.14e \n', param);
    fprintf(fid, '%.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f, %.1f \n', param);
    
end
    
% Beam Elastic Material
% fprintf(fid, '*Material, name=BeamMatl \n');
% fprintf(fid, '*Hyperelastic, neo hooke \n');
% fprintf(fid, '%2.3e, %2.3e \n', [2.5*mod_ratio, 0.6/mod_ratio]);
    
% Beam Elastic Material
fprintf(fid, '*Material, name=BeamMatl \n');
fprintf(fid, '*Elastic \n');
fprintf(fid, '%2.5e, %2.5e \n', param_beam);


% TIME POINTS
fprintf(fid,'*TIME POINTS, name=MustPoints, GENERATE\n');
fprintf(fid,'0., 1., %f\n',1/SolverOptions.n_incr);
%% STEP
InitStep = SolverOptions.InitStep;
minStep = SolverOptions.minStep;
maxStep = SolverOptions.maxStep;
StepDisp = linspace(0, X_disp,nSteps);

for i = 1:nSteps
    fprintf(fid, '*Step, INC=1000, name=Step-%d, nlgeom=YES \n',i);
    if SolverOptions.STBL
        fprintf(fid, '*Static, STABILIZE=%.4e\n',SolverOptions.STBLfac);
    else
        fprintf(fid, '*Static \n');
    end
   % fprintf(fid, '*Static, STABILIZE, factor=1e-5\n');
    
    fprintf(fid, '%.2e, 1., %.2e, %.2e \n', [InitStep, minStep, maxStep]);

    fprintf(fid, '*Boundary \n');
 
        
%         if strcmp(RPL, 'SSxy')
%         
%         fprintf(fid, 'FreeRefNode, 1,1,%.3e \n',StepDisp(i));
%         fprintf(fid, 'FreeRefNode, 2,6\n');
%         
%         elseif strcmp(RPL, 'SSyz')
%         
%         fprintf(fid, 'FreeRefNode, 1,2\n');
%         fprintf(fid, 'FreeRefNode, 3,3,%.3e \n',StepDisp(i));
%         fprintf(fid, 'FreeRefNode, 4,6\n');
%         end
%         
 
    if strcmp(RPL, 'UEx')
        fprintf(fid, '%d, 1, 1,%.3e \n',[ntotal+4,StepDisp(i)]);
    elseif strcmp(RPL, 'UEy')
        fprintf(fid, '%d, 2, 2,%.3e \n',[ntotal+5,StepDisp(i)]);
    elseif strcmp(RPL, 'UEz')
        fprintf(fid, '%d, 3, 3,%.3e \n',[ntotal+6,StepDisp(i)]);
    elseif strcmp(RPL, 'SSxy')
        fprintf(fid, '%d, 1, 1,%.3e \n',[ntotal+4,StepDisp(i)]);
        fprintf(fid, '%d, 1, 3, 0 \n', ntotal+5);
        fprintf(fid, '%d, 1, 3, 0 \n', ntotal+6);
    elseif strcmp(RPL, 'SSxz')
        fprintf(fid, '%d, 2, 2,%.3e \n',[ntotal+5,StepDisp(i)]);
        fprintf(fid, '%d, 1, 3, 0 \n', ntotal+4);
        fprintf(fid, '%d, 1, 3, 0 \n', ntotal+6);
    elseif strcmp(RPL, 'SSyz')
        fprintf(fid, '%d, 3, 3,%.3e \n',[ntotal+6,StepDisp(i)]);
        fprintf(fid, '%d, 1, 3, 0 \n', ntotal+4);
        fprintf(fid, '%d, 1, 3, 0 \n', ntotal+5);   
        
    else 
        error('Displacement Mode NOT SUPPRORTED YET :)')
    end
    
    
    % Output
    fprintf(fid, '*Restart, write, frequency=0 \n');
%     fprintf(fid, '*Output, field, variable=PRESELECT \n');
%     fprintf(fid, '*Output, history\n');
%     fprintf(fid, '*Node Output, NSET=FreeRefNode\n');
%     fprintf(fid, 'RF1,RF2,RF3,U1 \n');
    fprintf(fid, '*Output, history, variable=PRESELECT, time points=MustPoints \n');
    fprintf(fid, '*Output, field, variable=PRESELECT, time points=MustPoints \n');
    fprintf(fid, '*ELEMENT Output\n');
    fprintf(fid, 'S,SE,SF,SM,SK \n');
    fprintf(fid, '*EL PRINT, ELSET=BrickElems\n');
    fprintf(fid, 'DG\n');
    fprintf(fid, '*End Step \n');
end

     
%% FUNCTIONS
    function A_set_srt = nodeRegister(A_set,B_set,Nodes_mat)
        nl = length(A_set);
        if nl~=length(B_set)
            error('Sth is wrong with nodesets');
        end
        
        % Loop through nodesets.
        A_xyz = Nodes_mat(A_set,:); 
        B_xyz = Nodes_mat(B_set,:); 
        IDSX = knnsearch(A_xyz,B_xyz);
        A_set_srt = A_set(IDSX);
        
        % Check; If register is ok, all vectors should be equal, to the
        % difference between periodic faces/edges
        A_xyz_new = Nodes_mat(A_set_srt,:); 
        diff_vec  = A_xyz_new-B_xyz;
        checkIDs = ismember(diff_vec,diff_vec(1,:),'rows')';
        if sum(checkIDs)<nl
            error('Registering is worng')
        end
              
    

    end
        

    function []=writeNSET(fid_2, NSETx, NSET_name_x)

        n_in_set = length(NSETx);
        fprintf(fid_2, ['*NSET, NSET=',NSET_name_x,', unsorted \n']);
        if n_in_set < 8
            str_line_c = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
            fprintf(fid_2, str_line_c, NSETx);
        else
            n8rows_c = fix(n_in_set/8);
            nrest_c = mod(n_in_set,8);
            str_line_c = ['%i', repmat(', %i',1,7), ', \n'];
            fprintf(fid_2, str_line_c, NSETx(1:8*n8rows_c));
            if nrest_c>0
                str_line_c = ['%i', repmat(', %i',1,nrest_c-1), ' \n'];
                fprintf(fid_2, str_line_c, NSETx(8*n8rows_c+1:end));
            end
        end
    end

    function []=writeCONSTRAINT_UE(fid3,NSET_a, NSET_b, RP_i,A_coef,B_coef,C_coef,ntotal_sol,ASET,BSET)
        
        for ij = 1:3
            
         for ijj = 1:length(ASET)
                
            fprintf(fid3, ['%d,%d,%d\n'],[ASET(ijj),ij,A_coef]);
            fprintf(fid3, ['%d,%d,%d\n'],[BSET(ijj),ij,B_coef]);
            if C_coef(ij)~=0
                if length(RP_i)==1
                    fprintf(fid3, '%d,%d,%d\n',[ntotal_sol+RP_i,ij,C_coef(ij)]);
                else
                    fprintf(fid3, '%d,%d,%d\n',[ntotal_sol+RP_i(ij),ij,C_coef(ij)]);
                end
            end
         end
     
        end
    end



    function []=writeCONSTRAINT_SS(fid3,NSET_a, NSET_b, RP_i,RP_j,RP_k,i_coef,j_coef,k_coef,ntotal_sol,ASET,BSET)
        
        for ij = 1:3 
            for ijj = 1:length(ASET)
                fprintf(fid3, ['%d,%d,%d\n'],[ASET(ijj),ij,1]);
                fprintf(fid3, ['%d,%d,%d\n'],[BSET(ijj),ij,-1]);
                fprintf(fid3, '%d,%d,%d\n',[ntotal_sol+RP_i(ij), ij, i_coef]);
                fprintf(fid3, '%d,%d,%d\n',[ntotal_sol+RP_j(ij), ij, j_coef]);
                fprintf(fid3, '%d,%d,%d\n',[ntotal_sol+RP_k(ij), ij, k_coef]);
            end
        end
    end






end