%% Main Script for Beam-to-solid coupling modeling
% Based on the work by Steinbrecher et al., 2020
% Please follow tutorial on GitHub for more info. Thanks!
% Sotiris Kakaletsis, 2022-2023
clear; close all;
addpath('./Functions')

%% ======================== INPUT =========================================
%------- CONFIGURATION OPTIONS --------------------------------------------
JobNum = 9024; % Abaqus Job Number [0-9999]
JobDescription = 'Helix Coil Example'; % Analysis Description comment
WorkDir = './AbaqusWorkDir/'; % Local work directory for Abaqus
MatrixPaths = '/home/datastore/Sotiris/'; % Directory to save matrix *.txt

%------- GEOMETRY OPTIONS -------------------------------------------------
% Solid Domain Dimensions
W = 1;  % Global X-direction
D = 2;  % Global Y-direction
H = 1;  % Global Z-direction

% Beam/Fiber Geometry
r0 = [W/2, 0.05*D, H/2]; % Start point
rf = [W/2, 0.95*D, H/2]; % End point

% Discretization
hsolid = 1/18; % Solid element length
hbeam = 1.0*hsolid; % Beam element length
dist_tol = 1e-8; % Distance tolerance for merging beam nodes

% Introduce Undulations
curve_type='helix'; % 'straight' OR 'helix' OR 'sin_und'
Ax = min([W,H])*0.4; Ay=[]; % Undulations Amplitude
wx=[]; wy=[];  % Unudations frequency (applicable only for 'sin_und'
nloops=3; % Complete loops (only for 'helix')
BeamAxes = {r0, rf, Ax, Ay, wx, wy, nloops, curve_type, hbeam}; % Gather

%------- CONSTITUTIVE MODEL OPTIONS ---------------------------------------
PenaltyConst = 1e-2; % Coupling penalty parameter
nGP = 6; % Number of Gauss points for beam domain integration
 
% Material and geometry
BeamOrder = 2; % Linear OR quadratic Timoshenko Beams OR Euler-Bernoulli 
BeamRadius = 0.01;
mod_ratio = 5000; % Beam-to-solid stiffness ratio
Ebeam  = 6.5; % Young's modulus for beam material
beam_nu = 0.495; % Poisson's Ratio for beams


% Solid Material Parameters
Matl = 'NHK'; %  NHK (Neo Hookean) 
kmurat = 1e3; % Bulk-to-shear modulus ratio (nearly incompressible)

%------- BVP AND SOLVER OPTIONS -------------------------------------------
% Displacement Mode
RPL = 'SSxy'; % Displacement mode
X_disp = Lm; % Prescribed displacement

% Abaqus Solver Options Steps
nSteps = 100;
SolverOptions.InitStep = 1e-2;% Initial step
SolverOptions.minStep = 1e-4; % Min step
SolverOptions.maxStep = 1e-2; % Max step
SolverOptions.n_incr = 100; % Number of increments to ouput results
SolverOptions.STBL = true; % Stabilization boolean
SolverOptions.STBLfac = 2e-4; % If on, stabilization factor

%------- ADVANCED PIPELINE SETTINGS ---------------------------------------
ReDoDiscretization = true; % Boolean to re-calculate coupling matrices
LMOrder = BeamOrder; % Lagrange Multiplier Order, equal to beam order
s_resol = 6; % Segmentation resolution per beam element
COLim = 1e-7; % Cut off limit for  stiffness matrix entries
n_group = 3; % Group DOFs at super element
plot_segmentation  = false; % Boolean to inspect segmentation
nplotbeams = 7; % Beam domain interpolation data points for post-processing
LM_bool = true; % Plot LM field

%--------------------------------------------------------------------------
%---------------------- END OF INPUT --------------------------------------
%--------------------------------------------------------------------------
%% Beam Cross Section Data
BCS.A = pi*BeamRadius^2;
BCS.I = pi/4*BeamRadius^4;
BCS.J = pi/2*BeamRadius^4;
BCS.E = param_beam(1);
BCS.nu = param_beam(2);
BCS.mu = param_beam(1)/(2*(1+BCS.nu));

if strcmp(Matl, 'NHK')
    param_beam(1) = Ebeam;
    param_beam(2) = beam_nu;
        
    param_solid(1) = param_beam(1)/(4*mod_ratio*(1+beam_nu)); %C10
    param_solid(2) = 1/(param_solid(1)*kmurat);
    
end

%% DISCRETIZATION
if ReDoDiscretization
    
    % Beam Domain
    fprintf('STEP 1: Generating BEAM mesh...')
    [BeamMesh] = BeamMeshGenerator_v2(BeamAxes, BeamOrder);
    fprintf('DONE! \n')
     
    % Solid Domain
    nx = round(W/hsolid); 
    ny = round(D/hsolid); 
    nz = round(H/hsolid); 
    
    fprintf('STEP 2: Generating SOLID mesh...')
    SolidMesh = CubeMeshGenerator(nx,ny,nz, W, D, H);
    fprintf('DONE! \n')  

    % Find Common Beam to Solid Segments.
    fprintf('STEP 3: Beam-Solid SEGMENTATION...')
    [BeamMesh, npmax] = BeamSegmentGen_v3(W, D, H, nx, ny, nz, SolidMesh, BeamMesh, s_resol, plot_segmentation);
    fprintf('DONE! \n')

    % Create global connectivity matrix
    % Merge Connectivity and nodes matrix. Also, merge nodes (cross links). Define n1 and n2. 
    [GlobalBeamMesh] = GlobalBeamMeshGen_v2(BeamMesh, dist_tol);
    GlobalBeamMesh.npmax = npmax;
    GlobalBeamMesh.Order = BeamOrder;

    % EXPORT ALL DATA
    % Make an new working directory
    AbaqusNewDir = [WorkDir,'Job',num2str(JobNum,'%.4d'),'/'];
    if ~exist(AbaqusNewDir, 'dir')
       mkdir(AbaqusNewDir)
    end



    %% Mortar matrices
    fprintf('\nSTEP 4: Itegrating mortar matrices...\n')

    % Integrate and Assemble matrix M
    [M_g] = AssembleM_v3([], GlobalBeamMesh,SolidMesh, LMOrder, nGP);

    % Integrate D and k matrices
    [D_g, k_g_inv] = AssembleDnK_v3(GlobalBeamMesh, LMOrder, nGP);

    % Clear round of errors
    M_g(abs(M_g)<1e-20) = 0;
    D_g(abs(D_g)<1e-20) = 0;
    %k_g_inv(abs(k_g_inv)<1e-12) = 0;
    fprintf('DONE! \n')

    %% GROUP DOFs
    fprintf('STEP 5: Assembly SUPER ELEMENT...\n')
    MkM = (M_g'*(k_g_inv*M_g));
    MkD = (M_g'*(k_g_inv*D_g));
    DkM = MkD';
    %DkM = (D_g'*(k_g_inv*M_g));
    DkD = (D_g'*(k_g_inv*D_g));
    nS = size(SolidMesh.Nodes,1);
    nB = size(GlobalBeamMesh.Nodes,1);

    A_temp =[MkM,-MkD;-DkM,DkD];

    A_nodal = A_temp(1:3:end,1:3:end);
    cut_off_term = max(max(abs(A_nodal)))*COLim;
    A_nodal(abs(A_nodal)<cut_off_term) = 0; % Clean up the round off errors.

    % First, check solid nodes ds
    %Find Tempoarry connectivity
    for i = 1:nS+nB
        [~,J_nodes,~]=find(A_nodal(i,:));

        % Independent node should always be included in the node list
        if isempty(find(i==J_nodes)) && ~isempty(J_nodes)
            Jtemp = J_nodes;
            J_nodes = [Jtemp(Jtemp<i),i,Jtemp(Jtemp>i)];
            %disp('Relax round_off error correction')
        end

        % Add to collection of coupled nodes
        coupled_nodes{i,1} = J_nodes;

        n_cpl_nodes(i) = length(J_nodes);
        if n_cpl_nodes(i)>0

        end

    end

    
    % Sort coupled DOFs according to how many active nodes involve
    [~,IDs] = sort(n_cpl_nodes,'descend');
    IDs(n_cpl_nodes(IDs)==0)=[]; % Remove the inactive nodes
    DOFs_pool = IDs;

    % Initialize

    stop_grouping = false;
    c_grp = 0; % Group counter
    actv_nodes = {};
    n_actv_nodes = [];
    new_grps = [];
    while ~stop_grouping

        i_c = DOFs_pool(1); % Candidate DOF to be grouped
        DOFs_pool(1) = []; % Remove the initial candidate from the pool
        c_grp = c_grp +1;

        n_common = [];
        nodes2Bmatched = coupled_nodes{i_c,1};
        for i = 1:length(DOFs_pool)
            nodes2match = coupled_nodes{DOFs_pool(i),1};
            n_common(i) = length(intersect(nodes2Bmatched,nodes2match));
        end
        [~,IDs2choose] = sort(n_common,'descend');
        new_grps(c_grp,:) = [i_c,DOFs_pool(IDs2choose(1:n_group-1))];

        % Calculate active nodes in the group
        temp = [];
        for ii = 1:n_group
            temp = union(temp,coupled_nodes{new_grps(c_grp,ii),1});
        end
        actv_nodes{c_grp} = temp;
        n_actv_nodes(c_grp) = length(temp);

        % Erase chosen DOFs from the pool
        DOFs_pool(IDs2choose(1:n_group-1)) = [];

        % if none left
        if length(DOFs_pool)<=n_group-1
            stop_grouping = true;
        end

    disp('')
    end


    % Complete DOFs_pool
    if mod(length(IDs),n_group)~=0
        c_grp = c_grp + 1;
        new_grps(c_grp,:) = [DOFs_pool, zeros(1,n_group-length(DOFs_pool))];
        % Calculate active nodes in the group
        temp = [];
        for ii = 1:length(DOFs_pool)
            temp = union(temp,coupled_nodes{new_grps(c_grp,ii),1});
        end
        actv_nodes{c_grp} = temp;
        n_actv_nodes(c_grp) = length(temp);
    end

    fprintf('DONE! \n')
    %% GATHER, ARRANGE and COMPLETE MATRICES
    % Consider two groups
    nmax_AB = max(n_actv_nodes);
    A_Group = find(n_actv_nodes>1.5*nmax_AB/3);
    B_Group = find(n_actv_nodes<=1.5*nmax_AB/3);
    C_Group = find(n_actv_nodes<=nmax_AB/8 );

    BCcommon = intersect(B_Group,C_Group);
    B_Group = setxor(B_Group,BCcommon);

    figure()
    hold on
    plot(sort(n_actv_nodes,'descend'))
    yline(1.5*nmax_AB/3,'--')
    yline(nmax_AB/8,'--')
    hold off
    xlabel('DOFs')
    ylabel('Number of Active Nodes per DOF')

    %-------CONSIDER A GROUP-----------------------------------
    nmax_A = max(n_actv_nodes(A_Group));
    ncpl_A = length(A_Group); % number of coupling elements
    nT_A = ncpl_A*n_group;

    % Initialize
    CoupleElemCon_A = zeros(ncpl_A,nmax_A);
    AMTRX = zeros(nT_A,nmax_A);

    % Loop each coupling element
    for iel = 1:ncpl_A
        i = A_Group(iel);
        % Export DOFs
        DOFs_ORG = new_grps(i,:)';
        DOFs = new_grps(i,:)';

        % Export active nodes
        temp_NLIST = actv_nodes{i};

        % Re-arrange NLIST so that DOFs come first
        temp = ismember(temp_NLIST,DOFs);
        ID_dofs = find(temp);
        ID_rest = find(~temp);
        NLIST = [temp_NLIST(ID_dofs);temp_NLIST(ID_rest)]';
        DOFs = temp_NLIST(ID_dofs); % To ensure the order matches

        n_NLIST = length(NLIST); % Number of active nodes
        n_DUMMY = nmax_A-n_NLIST; % Nodes to be added, with 0 coefficient

        % Construct Local
        r_id = (iel-1)*n_group+1;
        if isempty(DOFs_ORG(DOFs_ORG==0))
            AMTRX(r_id:r_id+n_group-1,1:n_NLIST) = A_nodal(DOFs,NLIST);
        else
            % Last element, contains some zeros
            for ii = 1:length(DOFs)
                AMTRX(r_id+ii-1,1:n_NLIST) = A_nodal(DOFs(ii),NLIST);
            end
        end

        % Fill the Coupling Element connectivity with dummy nodes
        temp = setxor(NLIST,1:nS);
        rnd_ind = randperm(length(temp),n_DUMMY);
        CoupleElemCon_A(iel,:) = [NLIST,temp(rnd_ind)];          

    end


    %-------CONSIDER B GROUP----------------------------------
    nmax_B = max(n_actv_nodes(B_Group));
    ncpl_B = length(B_Group); % number of coupling elements
    nT_B = ncpl_B*n_group;

    % Initialize
    CoupleElemCon_B = zeros(ncpl_B,nmax_B);
    BMTRX = zeros(nT_B,nmax_B);

    % Loop each coupling element
    for iel = 1:ncpl_B
        i = B_Group(iel);
        % Export DOFs
        DOFs_ORG = new_grps(i,:)';
        DOFs = new_grps(i,:)';

        % Export active nodes
        temp_NLIST = actv_nodes{i};

        % Re-arrange NLIST so that DOFs come first
        temp = ismember(temp_NLIST,DOFs);
        ID_dofs = find(temp);
        ID_rest = find(~temp);
        NLIST = [temp_NLIST(ID_dofs);temp_NLIST(ID_rest)]';
        DOFs = temp_NLIST(ID_dofs); % To ensure the order matches

        n_NLIST = length(NLIST); % Number of active nodes
        n_DUMMY = nmax_B-n_NLIST; % Nodes to be added, with 0 coefficient

        % Construct Local
        r_id = (iel-1)*n_group+1;
        if isempty(DOFs_ORG(DOFs_ORG==0))
            BMTRX(r_id:r_id+n_group-1,1:n_NLIST) = A_nodal(DOFs,NLIST);
        else
            % Last element, contains some zeros
            for ii = 1:length(DOFs)
                BMTRX(r_id+ii-1,1:n_NLIST) = A_nodal(DOFs(ii),NLIST);
            end
        end

        % Fill the Coupling Element connectivity with dummy nodes
        temp = setxor(NLIST,1:nS);
        rnd_ind =  randperm(length(temp),n_DUMMY);
        CoupleElemCon_B(iel,:) = [NLIST,temp(rnd_ind)];

    end


    %-------CONSIDER C GROUP----------------------------------
    nmax_C = max(n_actv_nodes(C_Group));
    ncpl_C = length(C_Group); % number of coupling elements
    nT_C = ncpl_C*n_group;

    % Initialize
    CoupleElemCon_C = zeros(ncpl_C,nmax_C);
    CMTRX = zeros(nT_C,nmax_C);

    % Loop each coupling element
    for iel = 1:ncpl_C
        i = C_Group(iel);
        % Export DOFs
        DOFs_ORG = new_grps(i,:)';
        DOFs = new_grps(i,:)';

        % Export active nodes
        temp_NLIST = actv_nodes{i};

        % Re-arrange NLIST so that DOFs come first
        temp = ismember(temp_NLIST,DOFs);
        ID_dofs = find(temp);
        ID_rest = find(~temp);
        NLIST = [temp_NLIST(ID_dofs);temp_NLIST(ID_rest)]';
        DOFs = temp_NLIST(ID_dofs); % To ensure the order matches

        n_NLIST = length(NLIST); % Number of active nodes
        n_DUMMY = nmax_C-n_NLIST; % Nodes to be added, with 0 coefficient

        % Construct Local
        r_id = (iel-1)*n_group+1;
        if isempty(DOFs_ORG(DOFs_ORG==0))
            CMTRX(r_id:r_id+n_group-1,1:n_NLIST) = A_nodal(DOFs,NLIST);
        else
            % Last element, contains some zeros
            for ii = 1:length(DOFs)
                CMTRX(r_id+ii-1,1:n_NLIST) = A_nodal(DOFs(ii),NLIST);
            end
        end

        % Fill the Coupling Element connectivity with dummy nodes
        temp = setxor(NLIST,1:nS);
        rnd_ind =  randperm(length(temp),n_DUMMY);
        CoupleElemCon_C(iel,:) = [NLIST,temp(rnd_ind)];

    end


    %% Export Mortar matrices in abaqus form
    fprintf('STEP 6: Writing Mortar Matrices...\n')
    WriteFortranMatrix(AMTRX, MatrixPaths, 'A', JobNum)
    WriteFortranMatrix(BMTRX, MatrixPaths, 'B', JobNum)
    WriteFortranMatrix(CMTRX, MatrixPaths, 'C', JobNum)
    fprintf('DONE! \n')


    %% UPDATE UEL PARAMETERS
    fprintf('STEP 6: Updating UEL...')
    UpdateUEL_v6(WorkDir,nT_A,nT_B,nT_C,nmax_A,nmax_B,nmax_C,n_group, MatrixPaths, JobNum)
    fprintf('DONE! \n')


    CoupleElemCon.A = CoupleElemCon_A;
    CoupleElemCon.B = CoupleElemCon_B;
    CoupleElemCon.C = CoupleElemCon_C;
    nmax.A = nmax_A;
    nmax.B = nmax_B;
    nmax.C = nmax_C;
    nT.A = nT_A;
    nT.B = nT_B;
    nT.C = nT_C;        

    %Export current workspase with all data
    save([AbaqusNewDir,'SimData'])

else

    LocSimData = [WorkDir,'Job',num2str(JobNum,'%.4d'), '/SimData.mat'];
    load(LocSimData,'CoupleElemCon','nT','nmax','M_g','D_g','k_g_inv','SolidMesh', 'GlobalBeamMesh');

end

%% WRITE ABAQUS INPUT FILE
fprintf('STEP 7: Writing Abaqus input file...')
D01_WriteAbaqusInput_RPL_PBC(WorkDir,JobNum,SolidMesh, GlobalBeamMesh, PenaltyConst, BeamRadius, param_beam,  X_disp, nSteps, Matl, param_solid, beam_nu, CoupleElemCon, nmax, nT, SolverOptions,RPL)
fprintf('DONE! \n')
fprintf('----Matlab post processing COMPLETED--- \n')


%% RUN ABAQUS FILE
num_threads = 32;
cur_path = pwd;
IDStr = num2str(JobNum,'%.4d');

full_workdir = [cur_path,'/AbaqusWorkDir/Job',IDStr,'/'];
inp_file = [full_workdir, 'Model',num2str(JobNum,'%.4d'),RPL, '.inp'];
user_file = [full_workdir, 'UEL',num2str(JobNum,'%.4d'), '.f']; 

%Delete old result reports if any
S_res_file = [full_workdir, 'Model',num2str(JobNum,'%.4d'),RPL, '_Fields.csv'];
RF_res_file = [full_workdir, 'Model',num2str(JobNum,'%.4d'), RPL,'_History.csv'];

if isfile(S_res_file);  delete(S_res_file); end
if isfile(RF_res_file);  delete(RF_res_file); end


cmd_str1 = ['cd ',full_workdir];
cmd_str2 = ['abaqus job=Model',IDStr,RPL,' user=UEL',IDStr,'.f double cpus=' num2str(num_threads), ' int > NUL ask_delete=OFF '];
%cmd_str3 = ['abaqus odbreport mode=csv odb=Model',IDStr,' job=Model',IDStr,'_S step=_LAST_ frame=_LAST_ field=S > NUL'];
frame_string = sprintf('%d,',[0:1:SolverOptions.n_incr+1]);
frame_string = frame_string(1:end-1);
cmd_str4 = ['abaqus odbreport mode=csv odb=Model',IDStr,RPL,' job=Model',IDStr,RPL,'_Fields  field=U,S,SE,SF,SM,SK,RF frame=',frame_string,' > NUL int'];
cmd_str5 = ['abaqus odbreport mode=csv odb=Model',IDStr,RPL,' job=Model',IDStr,RPL,'_History  history=ALLSE,ALLWK,ALLSD,ALLIE  > NUL int'];


tic
fprintf('STEP 8: ABAQUS running.')
[statusRun, blabla] = system([cmd_str1,'; ',cmd_str2]);
if statusRun ~= 0
    fprintf('Abaqus simulation failed.')
end
fprintf('...COMPLETE(Runtime %.1fs).\n\n\n',toc)

tic
fprintf('STEP 9: Generating ABAQUS reports.')
[statusRun, blabla] = system([cmd_str1,'; ',cmd_str4]);
[statusRun, blabla] = system([cmd_str1,'; ',cmd_str5]);
fprintf('...COMPLETE(Runtime %.1fs).\n\n\n',toc)


%% Create VTK files
fprintf('STEP 10: Importing ABAQUS reports...')
AbaqusData = ReadAbaqusReportsRPL(cur_path, JobNum,RPL);
fprintf('DONE! \n')
fprintf('STEP 11: Writing VTKs')
AbaqusData.NRGS = create_vtkRPL(cur_path, JobNum, AbaqusData, SolidMesh, GlobalBeamMesh, nplotbeams, PenaltyConst, M_g, D_g, k_g_inv,LM_bool,BCS,RPL);

%Save Rigid Force (Reaction)
RFx = AbaqusData.RFdisp;
RFy = AbaqusData.RF;
save([full_workdir,RPL,'AbaqusData.mat'],'AbaqusData')
fprintf('DONE! \n')


