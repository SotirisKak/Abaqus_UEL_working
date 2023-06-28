function [] = D08_WriteAbaqusInput_v6(WorkDir, JobNum, SolidMesh, disp_mode ,X_disp, nSteps, Matl, param_solid, SolverOptions)


filename = [WorkDir, '/Job',num2str(JobNum,'%.4d'),'/Model',num2str(JobNum,'%.4d'),'.inp'];
fid = fopen(filename,'w');

%% GEOMETRY SECTION
nS = size(SolidMesh.Nodes,1);

Nodes = [SolidMesh.Nodes];
ntotal = nS;
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
fprintf(fid, '\t %i, %2.14e, %2.14e, %2.14e \n', [ntotal+1, W/2, -0.1*D, H/2]);
fprintf(fid, '\t %i, %2.14e, %2.14e, %2.14e \n', [ntotal+2, W/2, 1.1*D, H/2]);
%% Solid Node Sets
FixSolid = SolidMesh.negx;
n_in_set = length(FixSolid);
fprintf(fid, '*NSET, NSET=NegX \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, FixSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, FixSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, FixSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.posx;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=PosX \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.posy;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=PosY \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.negy;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=NegY \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.posz;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=PosZ \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end

DispSolid = SolidMesh.negz;
n_in_set = length(DispSolid);
fprintf(fid, '*NSET, NSET=NegZ \n');
if n_in_set < 8
    str_line = ['%i', repmat(', %i',1,n_in_set-1), ' \n'];
    fprintf(fid, str_line, DispSolid);
else
    n8rows = fix(n_in_set/8);
    nrest = mod(n_in_set,8);
    str_line = ['%i', repmat(', %i',1,7), ', \n'];
    fprintf(fid, str_line, DispSolid(1:8*n8rows));
    if nrest>0
        str_line = ['%i', repmat(', %i',1,nrest-1), ' \n'];
        fprintf(fid, str_line, DispSolid(8*n8rows+1:end));
    end
end


%% Ref Node Sets
fprintf(fid, '*NSET, NSET=FixRefNode\n');
fprintf(fid, '%i \n', ntotal+1);
fprintf(fid, '*NSET, NSET=FreeRefNode\n');
fprintf(fid, '%i \n', ntotal+2);

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Xp_Face_S4\n');
fprintf(fid, '%i,\n', S4FaceBricks');
fprintf(fid, '*SURFACE, Name=Xp_Face\n');
fprintf(fid, 'Xp_Face_S4, S4\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Xn_Face_S6\n');
fprintf(fid, '%i,\n', S6FaceBricks');
fprintf(fid, '*SURFACE, Name=Xn_Face\n');
fprintf(fid, 'Xn_Face_S6, S6\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Yp_Face_S5\n');
fprintf(fid, '%i,\n', S5FaceBricks');
fprintf(fid, '*SURFACE, Name=Yp_Face\n');
fprintf(fid, 'Yp_Face_S5, S5\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Yn_Face_S3\n');
fprintf(fid, '%i,\n', S3FaceBricks');
fprintf(fid, '*SURFACE, Name=Yn_Face\n');
fprintf(fid, 'Yn_Face_S3, S3\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Zp_Face_S2\n');
fprintf(fid, '%i,\n', S2FaceBricks');
fprintf(fid, '*SURFACE, Name=Zp_Face\n');
fprintf(fid, 'Zp_Face_S2, S2\n');

% Element Set to define surface
fprintf(fid, '*ELSET, ELSET=Zn_Face_S1\n');
fprintf(fid, '%i,\n', S1FaceBricks');
fprintf(fid, '*SURFACE, Name=Zn_Face\n');
fprintf(fid, 'Zn_Face_S1, S1\n');

%% ElEMENT SECTION



% Bricks
fprintf(fid, '*ELEMENT, TYPE=C3D8RH, ELSET=BrickElems \n');
SCon = SolidMesh.Connectivity;
nSolEl = size(SCon,1); 
SolidConnectivity = [(1:nSolEl)',SCon];
fprintf(fid, ['%i', repmat(', %i',1,8), ' \n'], SolidConnectivity');



%% SECTIONS
fprintf(fid, '** \n');
fprintf(fid, '** SECTIONS \n');
fprintf(fid, '** \n');
if strcmp(Matl, 'NHK') || strcmp(Matl, 'OGD') || strcmp(Matl, 'MNRIV')
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


%% Constraints
fprintf(fid, '** \n');
fprintf(fid, '** Constraints \n');
fprintf(fid, '** \n');
fprintf(fid,'*Coupling, constraint name=FixedConstraint, ref node=FixRefNode, surface=Yn_Face\n');
fprintf(fid,'*Kinematic\n');
fprintf(fid,'*Coupling, constraint name=FreeConstraint, ref node=FreeRefNode, surface=Yp_Face\n');
fprintf(fid,'*Kinematic\n');
%% Materials
fprintf(fid, '** \n');
fprintf(fid, '** Materials \n');
fprintf(fid, '** \n');

if strcmp(Matl, 'NHK')%Solid Neo Hookean
    
    fprintf(fid, '*Material, name=SolidMatl \n');
    fprintf(fid, '*Hyperelastic, neo hooke \n');
    fprintf(fid, '%2.3e, %2.3e \n', param_solid);
    
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
    % Boundary Conditions
    fprintf(fid, '*Boundary \n');
    fprintf(fid, 'FixRefNode, ENCASTRE \n');
    fprintf(fid, '*Boundary \n');
    if strcmp(disp_mode, 'SimpleShear')
        
        fprintf(fid, 'FreeRefNode, 1,1,%.3e \n',StepDisp(i));
        fprintf(fid, 'FreeRefNode, 2,6\n');
        
    elseif strcmp(disp_mode, 'UniaxialExt')
        
        fprintf(fid, 'FreeRefNode, 1,1\n');
        fprintf(fid, 'FreeRefNode, 2,2,%.3e \n',StepDisp(i));
        fprintf(fid, 'FreeRefNode, 3,6\n');
        
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
    fprintf(fid, '*End Step \n');
    
end