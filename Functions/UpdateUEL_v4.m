function [] = UpdateUEL_v4(UELDir,nT,nGr,nmax, matpath, jobnum)

% Read txt into cell A
DefaultUELDir = [UELDir,'UEL_BeamCoupling_v4.f'];
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

A{23} = sprintf(['	      ADir = "',Afilename,'"']);


% Change values
A{15} = sprintf('	  PARAMETER (nT=%d,nG=%d,nM=%d)', [nT,nGr,nmax]);
A{46} = sprintf('	  PARAMETER (nT=%d,nG=%d,nM=%d)', [nT,nGr,nmax]);

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