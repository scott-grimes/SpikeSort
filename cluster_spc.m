function [clu, tree] = run_cluster(handles)
dim=handles.par.inputs;
fname=handles.par.fname;
fname_in=handles.par.fname_in;

% DELETE PREVIOUS FILES
fileexist = exist([fname '.dg_01.lab'],'file');
if(fileexist~=0)
    delete([fname '.dg_01.lab']);
    delete([fname '.dg_01']);
end

dat=load(fname_in);
n=length(dat);
fid=fopen(sprintf('%s.run',fname),'wt');
fprintf(fid,'NumberOfPoints: %s\n',num2str(n));
fprintf(fid,'DataFile: %s\n',fname_in);
fprintf(fid,'OutFile: %s\n',fname);
fprintf(fid,'Dimensions: %s\n',num2str(dim));
fprintf(fid,'MinTemp: %s\n',num2str(handles.par.mintemp));
fprintf(fid,'MaxTemp: %s\n',num2str(handles.par.maxtemp));
fprintf(fid,'TempStep: %s\n',num2str(handles.par.tempstep));
fprintf(fid,'SWCycles: %s\n',num2str(handles.par.SWCycles));
fprintf(fid,'KNearestNeighbours: %s\n',num2str(handles.par.KNearNeighb));
fprintf(fid,'MSTree|\n');
fprintf(fid,'DirectedGrowth|\n');
fprintf(fid,'SaveSuscept|\n');
fprintf(fid,'WriteLables|\n');
fprintf(fid,'WriteCorFile~\n');
if handles.par.randomseed ~= 0
    fprintf(fid,'ForceRandomSeed: %s\n',num2str(handles.par.randomseed));
end    
fclose(fid);

[str,maxsize,endian]=computer;
handles.par.system=str;
switch handles.par.system
    case {'PCWIN','PCWIN64'}    
        if exist([pwd '\cluster.exe'])==0
            directory = which('cluster.exe');
            copyfile(directory,pwd);
        end
        dos(sprintf('cluster.exe %s.run',fname));
end
        
clu=load([fname '.dg_01.lab']);
tree=load([fname '.dg_01']); 
delete(sprintf('%s.run',fname));
delete *.mag
delete *.edges
delete *.param
delete(fname_in); 
