function [channel_file_names] = extract_data(myfilename)
% Scott Grimes - Max Planck Cybernetics - 2011
% Extracts channels from a given .mat file for later use


newData1 = load(myfilename); %loads file for name extraction
vars = fieldnames(newData1); %gets channel names
%   for i = 1:length(vars)
%      assignin('caller', vars{i}, newData1.(vars{i}));
%   end
clear newData1
load(myfilename) %loads file for data extraction
for i = 1:length(vars)
fprintf('Extracting channel %i of %i...\n',i,length(vars));
myvarname = vars(i);
d = '.values';
L = '.length';
S = '.start';
I = '.interval';
d = strcat(myvarname,d);
L = strcat(myvarname,L);
S = strcat(myvarname,S);
I = strcat(myvarname,I);
try
eval(['data =' char(d) ';']); %loads data from channel=
eval(['L =' char(L) ';']); 
eval(['S =' char(S) ';']);
eval(['interval =' char(I) ';']); %loads interval from channel
eval('ts = [S:interval:interval*L]'';'); %builds time vector from channel
eval('data = data'';');
clear d L S I
catch ME
    channel_file_names = -100.1;
    return;
end
savefile = char(strcat('Data_',myvarname));
save(savefile,'data','interval','ts');

    
end
channel_file_names = strcat('Data_',vars')';