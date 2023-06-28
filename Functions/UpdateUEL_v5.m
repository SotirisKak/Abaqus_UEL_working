function [] = UpdateUEL_v5(UELDir,nT_A,nT_B,nmax_A,nmax_B,n_Gr, matpath, jobnum)

% Read txt into cell A
DefaultUELDir = [UELDir,'UEL_BeamCoupling_v5.f'];
fid = fopen(DefaultUELDir,'r');
i = 1;
tline = fgetl(fid);
A{i,1} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i,1} = tline;
end
fclose(fid);

% Replace matrix filenames
Afilename = [matpath,'A',num2str(jobnum,'%.4d'),'.txt'];
Bfilename = [matpath,'B',num2str(jobnum,'%.4d'),'.txt'];

A{23} = sprintf(['	ADir = "',Afilename,'"']);
A{24} = sprintf(['	BDir = "',Bfilename,'"']);



% Change values
A{15} = sprintf('      PARAMETER (nTA=%d,nMA=%d,nTB=%d,nMB=%d,nG=%d)', [nT_A,nmax_A,nT_B,nmax_B,n_Gr]);
A{48} = sprintf('      PARAMETER (nTA=%d,nMA=%d,nTB=%d,nMB=%d,nG=%d)', [nT_A,nmax_A,nT_B,nmax_B,n_Gr]);

AbaqusNewDir = [UELDir,'Job',num2str(jobnum,'%.4d')];
if ~exist(AbaqusNewDir, 'dir')
    mkdir(AbaqusNewDir)
end
UEL_NewDir = [AbaqusNewDir,'/UEL',num2str(jobnum,'%.4d'),'.f'];

% Write the new subroutine
fid = fopen(UEL_NewDir, 'w');
for i = 1:numel(A)
    if A{i+1} == -1
        fprintf(fid,'%s', A{i});
        break
    else
        fprintf(fid,'%s\n', A{i});
    end
end
fclose(fid);