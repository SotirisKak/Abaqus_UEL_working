function [] = create_vtk_nobeams(cur_path, JobNum, AbaqusData)

% Create paths needed
% Define paths
IDStr = num2str(JobNum,'%.4d');
full_workdir = [cur_path,'/AbaqusWorkDir/Job',IDStr,'/'];
inp_dir = [full_workdir,'Model',IDStr,'.inp'];
vtkDir = [full_workdir,'VTKs/'];
if ~exist(vtkDir, 'dir')
    mkdir(vtkDir)
else
    % Remove old results
    eval(['rmdir ' vtkDir ' s']);
    % Re-make directory
    mkdir(vtkDir)
end

solid_file = [vtkDir,'S_',IDStr];
beam_file = [vtkDir,'B_',IDStr];
seg_file = [vtkDir,'Seg_',IDStr];
elem_file = [vtkDir,'Elem_',IDStr];


%% Read Abaqus input file
fid2 = fopen(inp_dir,'r');
line2 = fgetl(fid2);

flagnode=0;
flagbrick=0;
flagtruss=0;
C3D8R = [];
C3D8I = [];
P = [];
B32 = [];
B31 = [];

while line2>=0

    if contains(line2,'*NODE')
        temp = fscanf(fid2,'%d,%f,%f,%f\n'); 
        temp = reshape(temp,[4 length(temp)/4 ]);
        temp = temp';
        P = [P;temp];
        temp = []; 
    end     
    
    if contains(line2,'C3D8R')
        temp = fscanf(fid2,'%d,%f,%f,%f,%f,%f,%f,%f,%f\n');
        temp = reshape(temp,[9 length(temp)/9]);
        temp = temp';
        C3D8R = [C3D8R; temp];
        temp = [];         
    end 
    
    if contains(line2,'C3D8I')
        temp = fscanf(fid2,'%d,%f,%f,%f,%f,%f,%f,%f,%f\n');
        temp = reshape(temp,[9 length(temp)/9]);
        temp = temp';
        C3D8I = [C3D8I; temp];
        temp = [];         
    end 
    if contains(line2,'B31')
        temp = fscanf(fid2,'%d,%f,%f\n');
        temp = reshape(temp,[3 length(temp)/3]);
        temp = temp';
        B31 = [B31;temp];
        temp = [];
    end 
    
    if contains(line2,'B32')
        temp = fscanf(fid2,'%d,%f,%f,%f\n');
        temp = reshape(temp,[4 length(temp)/4]);
        temp = temp';
        B32 = [B32;temp];
        temp = [];
    end 
    
    line2 = fgetl(fid2);
    
    if isempty(line2)
        line2=fgetl(fid2);
    end 

end 

%B32 = B32(~all(B32 == 0, 2),:); % Removes (0,0) entries from T3D2 

temp = C3D8R; 
for i=1:size(temp,1)
    C3D8R(temp(i,1),:) = temp(i,:);
end 

nS = length(unique(C3D8R(:,2:end)));


%% WRITE VTK Steps

n_steps = size(AbaqusData.U,3);
count = -1;

%Loop steps
for i = 1:n_steps
    
        count = count+1;
        Utemp=[];
        Utemp(:,:) = AbaqusData.U(:,:,i);
        Stemp(:,:) = AbaqusData.S(:,:,i);
        if max(Stemp(:,2))>1
            error('Solid Stress valid only for Reduced Integration Elements');
        end
        
        dS = Utemp(1:nS,2:4);
         %Flatten the global nodal values vectors
        dS = reshape(dS',3*size(dS,1),1);
       
        
        
        % Solid stress values
        sqmat = [];
        for ii = 1:size(Stemp,1)
            sqmat = [sqmat;diag(Stemp(ii,3:5)) + squareform(Stemp(ii,6:8))]; 
        end
                       
        
        %% Solid Elements
        file_ext = ['_',num2str(count,'%.5d'),'.vtk'];
        fid3 = fopen([solid_file,file_ext],'w');
        fprintf(fid3,'# vtk DataFile Version 2.0\n');
        fprintf(fid3,'test\n');
        fprintf(fid3,'ASCII\n');
        fprintf(fid3,'DATASET UNSTRUCTURED_GRID\n');

        fprintf(fid3,['POINTS ',num2str(size(P,1)),' float\n']);
        fprintf(fid3,'%f %f %f\n',(P(:,2:4)+ Utemp(:,2:4))');

        fprintf(fid3,['CELLS ',num2str(size(C3D8R,1)),' ',num2str(9*size(C3D8R,1)),'\n']);
        fprintf(fid3,'8 %d %d %d %d %d %d %d %d\n',C3D8R(:,2:9)'-1);

        fprintf(fid3,['CELL_TYPES ',num2str(size(C3D8R,1)),'\n']);
        fprintf(fid3,'%d\n',12*ones(1,size(C3D8R,1))); 
        
        fprintf(fid3,['CELL_DATA ',num2str(size(C3D8R,1)),'\n']);
        fprintf(fid3,'TENSORS S float\n');
        fprintf(fid3,'%f %f %f\n',sqmat');
        
         
        fprintf(fid3,['POINT_DATA ',num2str(size(P,1)),'\n']);
        fprintf(fid3,'VECTORS U float\n');
        fprintf(fid3,'%f %f %f\n',Utemp(:,2:4)');


        fclose(fid3);

         
        
    end
end
