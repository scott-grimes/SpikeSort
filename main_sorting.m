function varargout = main_sorting(varargin)
% Scott Grimes - Max Planck Cybernetics - 2011
% Spike sorting software
% MAIN_SORTING M-file for main_sorting.fig
%      MAIN_SORTING, by itself, creates a new MAIN_SORTING or raises the existing
%      singleton*.
%
%      H = MAIN_SORTING returns the handle to a new MAIN_SORTING or the handle to
%      the existing singleton*.
%
%      MAIN_SORTING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_SORTING.M with the given input arguments.
%
%      MAIN_SORTING('Property','Value',...) creates a new MAIN_SORTING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_sorting_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_sorting_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_sorting

% Last Modified by GUIDE v2.5 14-Jul-2011 11:54:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_sorting_OpeningFcn, ...
                   'gui_OutputFcn',  @main_sorting_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main_sorting is made visible.
function main_sorting_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_sorting (see VARARGIN)

% Choose default command line output for main_sorting
cyberlogo = imread('logo.png');
axes(handles.axes8);
image(cyberlogo);
axis off
set(handles.checkbox_for_plotting_all,'Value',1);


handles.output = hObject;

% Update handles structure
set(hObject,'toolbar','figure');
guidata(hObject, handles);

% UIWAIT makes main_sorting wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_sorting_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
guidata(hObject, handles); 


function main_file_name_text_Callback(hObject, eventdata, handles)
% hObject    handle to main_file_name_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of main_file_name_text as text
%        str2double(get(hObject,'String')) returns contents of main_file_name_text as a double
guidata(hObject, handles); 

% --- Executes during object creation, after setting all properties.
function main_file_name_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to main_file_name_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.



if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); 

% --- Executes on button press in button_main_file_name. %file name to be loaded
function button_main_file_name_Callback(hObject, eventdata, handles)
% hObject    handle to button_main_file_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
logf = fopen('logfile.txt','w');fclose(logf); diary logfile.txt;
myfilename = get(handles.main_file_name_text,'String'); %gets filename to be used
channel_file_names = extract_data(myfilename); %gets channels from filename
% if iscell(channel_file_names)==1
%     for i=1:length(channel_file_names)
%     channel_file_name1(i,:) = channel_file_names{i};
%     end
% channel_file_names = channel_file_name1;
% end
save('deleteme.mat','channel_file_names')
if iscell(channel_file_names) ~= 1
    e = 1;
    save('invalidload.mat','e')
    fprintf('File not exported from spike2 correctly, breaking...\n')
else
cutoff = get(handles.cutoff_textbox,'String'); %gets cutoff
cutoff = str2num(cutoff);
milisecleft = get(handles.spike_length_miliseconds_left,'String'); %gets spike length in ms for left side
milisecleft = str2num(milisecleft);
milisecright = get(handles.spikelengthmiliseconds_right,'String'); %gets spike length in ms for right side
milisecright = str2num(milisecright);

for i = 1:length(channel_file_names(:,1))
    
fprintf('Processing channel %i of %i...\n',i,length(channel_file_names(:,1)))
try

channel_being_loaded = channel_file_names{i,:};
catch del1
channel_being_loaded = channel_file_names(i,:);
end
load(channel_being_loaded); %loads each channel


[spikes_temp spike_peak_time_temp] = obtain_spikes(data,interval,ts,0,cutoff,milisecleft,milisecright); %obtains spikes for each channel
if i ==1
    spikes = spikes_temp; %spikes from each channel added to eachother
    spike_peak_time = spike_peak_time_temp; %spike timestamps added to each other
else
    spikes = [spikes; spikes_temp]; %spikes from each channel added to eachother
    spike_peak_time = [spike_peak_time; spike_peak_time_temp]; %spike timestamps added to each other
end
clear spikes_temp %clears current channels spikes
delete(strcat(char(channel_being_loaded),'.mat')); %removes channel after spikes have been loaded
end
load('ansc1.mat')
fprintf('Aligning Spikes...\n');
clear data filtered



%removing nullspikes
nulls = find(sum(spikes')~=0); %finds null spikes
spikes = spikes(nulls,:); %removes nullspikes
spike_peak_time = spike_peak_time(nulls,:); %removes timestamps of nullspikes
spike_num = length(spikes(:,1)); %new number of spikes
spikes = spikes(:,find(sum(abs(spikes))~=0)); %removes nulls at begining/end
clear nulls
n = zeros(spike_num,1);
%save('tempfile.mat','spikes')
%finds max value cell index 
  for i=1:spike_num;
      n(i,1)= find(spikes(i,:)==max(spikes(i,:)),1); %max value of every spike
  end
save('tstemp.mat','ts')
clear ts
 %centers signals based on max value
 center=max(n);
 spikes_shifted=[zeros(1,length(spikes(1,:))*2)];
 for i=1:length(spikes(:,1));
     shift_num=center-n(i,1);
     if shift_num>0
         %shifts spikes so the abs(max) value of each spike is in the same
         %column for all spikes
     spikes_shifted(i,1:length(spikes_shifted(1,:)))=[zeros(1,shift_num),spikes(i,:),zeros(1,length(spikes_shifted(1,:))-shift_num-length(spikes(i,:)))];
     end
 end
fprintf('100%% Complete\n%i Spikes Found\n',spike_num)
%truncates null signals before/after spike
 spikes_shifted=spikes_shifted(:,find((max(spikes_shifted)==0)==0)); %cleans up the beginning and end of the spikes
 spikes=(spikes_shifted);
 clear spikes_shifted  n shift_num center

load tstemp
delete sorted_spikes.mat tstemp.mat
savefile = 'sorted_spikes.mat';
idx = zeros(size(spike_peak_time));
save(savefile,'ts','spikes','spike_peak_time','idx','interval')
savefile2 = 'ansc1.mat';
load(savefile2,'filtered')
save(savefile2,'spike_num', 'spike_peak_time','cutoff','filtered','s','m','amp_thresh','interval')




fprintf('Spikes Saved!\n');
clear data spikes 
end
guidata(hObject, handles); 

function channels_to_process_Callback(hObject, eventdata, handles)
% hObject    handle to channels_to_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channels_to_process as text
%        str2double(get(hObject,'String')) returns contents of channels_to_process as a double


% --- Executes during object creation, after setting all properties.
function channels_to_process_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channels_to_process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_process_channels.
function button_process_channels_Callback(hObject, eventdata, handles)
% hObject    handle to button_process_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles); 



% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_to_initiate_clustering.
function button_to_initiate_clustering_Callback(hObject, eventdata, handles)
% hObject    handle to button_to_initiate_clustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(handles.listbox2,'Value')
    case 1
        %using PCA clustering
   load sorted_spikes.mat
optional_cluster_override_number = get(handles.text_box_for_optional_cluster_override,'String');
optional_cluster_override_number = str2num(optional_cluster_override_number);
idx=cluster(spikes,optional_cluster_override_number,1); % 0 for automatic cluster identification, or specify your own cluster number
k = max(idx);
justthefacts = 'sorted_spikes.mat';
save(justthefacts,'spikes','idx','ts','spike_peak_time','interval');
clear justthefacts ans
    case 2
        %using SPC clustering
        load sorted_spikes.mat
        index = zeros(1,length(spikes(:,1)));
        save('tempSPCfile.mat','spikes','ts','index')
fout=fopen('files.txt', 'w+');
    fprintf(fout,'tempSPCfile.mat');
    fclose(fout);
        figure
    do_clustering()
    load times_tempSPCfile.mat
    idx = cluster_class(:,1);
clc
fprintf('SPC Clustering Complete: Found %i Clusters\n',max(idx));
    justthefacts = 'sorted_spikes.mat';
save(justthefacts,'spikes','idx','ts','spike_peak_time','interval');
clear justthefacts ans cluster_class index inspk ipermut par spikes ts
delete times_tempSPCfile.mat tempSPCfile.mat data_tempSPCfile.mat.dg_01 data_tempSPCfile.mat.dg_01.lab cluster_results.txt files.txt
    case 3 
    %using custom clustering
load sorted_spikes.mat
optional_cluster_override_number = get(handles.text_box_for_optional_cluster_override,'String'); %gets optinal cluster number
optional_cluster_override_number = str2num(optional_cluster_override_number);
idx=cluster(spikes,optional_cluster_override_number,0); % 0 for automatic cluster identification, or specify your own cluster number
k = max(idx); %number of clusters
justthefacts = 'sorted_spikes.mat';
save(justthefacts,'spikes','idx','ts','spike_peak_time','interval');
clear justthefacts ans
    
        
        
    otherwise
end
fprintf('Clustering Complete!\n');
set(handles.clusters_to_plot_textbox,'String','0');
set(handles.checkbox_for_plotting_all,'Value',1);
guidata(hObject, handles); %updates the handles
plot_seperated_clusters_button_Callback(hObject, eventdata, handles)
guidata(hObject, handles); 


function clusters_to_plot_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to clusters_to_plot_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of clusters_to_plot_textbox as text
%        str2double(get(hObject,'String')) returns contents of clusters_to_plot_textbox as a double
guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function clusters_to_plot_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to clusters_to_plot_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); %updates the handles

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plot_seperated_clusters_button.
function plot_seperated_clusters_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_seperated_clusters_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load sorted_spikes.mat
load('sorted_spikes.mat','interval')
colormatrix = colormap('Lines');

%axes(handles.axes_for_raw_data)

individual_spikes=figure;
load('sorted_spikes.mat','interval')
cla reset
hold on

clusters_to_plot = sscanf(get(handles.clusters_to_plot_textbox,'String'),'%i')'; %which clusters to be plot?
fastplot = get(handles.checkbox_for_plotting_all,'Value'); %plots only mean+std of clusters, use when clusters are very large to save time
if fastplot ==0
kk = zeros(15,3);
else
    kk = colormatrix;
end
sort(clusters_to_plot);
if clusters_to_plot ~=0
    tempidx = zeros(size(idx));
for i = 1:length(clusters_to_plot)
    z = clusters_to_plot(i);
    tempidx = tempidx+(idx==z);
end
spikes = spikes(find(tempidx==1),:);
idx = idx(find(tempidx==1));
end
clear tempidx z
if clusters_to_plot ~= 0
size_of_subplot = length(clusters_to_plot);
else
    size_of_subplot = max(idx);
    clusters_to_plot =1:max(idx);
end
%PLOTS SPIKES
k = max(idx);
if fastplot ==0
for i=1:length(spikes(:,1));
    subplot(size_of_subplot,1,find(clusters_to_plot==idx(i,1)))
    hold on
    plot(spikes(i,:),'Color',colormatrix(idx(i,1),:))
end
end

for kt = clusters_to_plot
    tempidx=kt==idx;
    temp = spikes(tempidx,:);
    temp = mean(temp);
    tempstd = std(spikes(tempidx,:));
    subplot(size_of_subplot,1,find(clusters_to_plot==kt))
    hold on
plot(temp,'Color',kk(kt,:),'LineWidth',2)
plot(temp+tempstd,'Color',kk(kt,:))
plot(temp-tempstd,'Color',kk(kt,:))
title(['Cluster ',num2str(kt)])
xlabel(['Points (',num2str(interval),'sec)'])
ylabel('mv')
clear temp tempstd
end
clear tempidx temp kt
%counts number in each cluster/titles for graphs
spike_num = length(spikes(:,1));
for i =clusters_to_plot
    fprintf('Cluster %i: %i Spikes (%%%i)\n',i,sum(idx==i),round(sum(idx==i)/spike_num*100))
end
fprintf('\n');
axes(handles.all_spikes_together_axis)
cla reset
title('Combined Spikes');
xlabel(['Points (',num2str(interval),'sec)'])
ylabel('mv')
hold on
if fastplot == 0
for i = 1:length(spikes(:,1));
    plot(spikes(i,:),'Color',colormatrix(idx(i,1),:))
end
else
    for kt = clusters_to_plot
    tempidx=kt==idx;
    temp = spikes(tempidx,:);
    temp = mean(temp);
    tempstd = std(spikes(tempidx,:));
plot(temp,'Color',kk(kt,:),'LineWidth',2)
plot(temp+tempstd,'Color',kk(kt,:))
plot(temp-tempstd,'Color',kk(kt,:))
clear temp tempstd
end
end
axes(handles.axes_for_just_cluster_averages)
cla reset
xlabel(['Points (',num2str(interval),'sec)'])
ylabel('mv')
hold on
for kt = 1:max(idx)
    tempidx=kt==idx;
    temp = spikes(tempidx,:);
    temp = mean(temp);
plot(temp,'Color',colormatrix(kt,:),'LineWidth',2)
clear temp
end

clear tempidx temp kt ts spikes spike_num size_of_subplot k idx colormatrix clusters_to_plot i individual_spikes ans




%plotting raw/filtered data




guidata(hObject, handles); %updates the handles


% --- Executes on selection change in timeline_menu_for_plotting.
function timeline_menu_for_plotting_Callback(hObject, eventdata, handles)
% hObject    handle to timeline_menu_for_plotting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns timeline_menu_for_plotting contents as cell array
%        contents{get(hObject,'Value')} returns selected item from timeline_menu_for_plotting


% --- Executes during object creation, after setting all properties.
function timeline_menu_for_plotting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeline_menu_for_plotting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_timeline_button.
function plot_timeline_button_Callback(hObject, eventdata, handles)
% hObject    handle to plot_timeline_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure
hold on
switch get(handles.timeline_menu_for_plotting,'Value')   
    case 1
        %plots the raw data
load('ansc1.mat','filtered','cutoff')
plot(filtered)
title('Data')
m = mean(filtered);
s = std(abs(filtered));
plot([0,length(filtered)],[m+cutoff*s,m+cutoff*s],'-r')
plot([0,length(filtered)],[m-cutoff*s,m-cutoff*s],'-r')
clear filtered s m amp_thresh spike_num ts spike_peak_time spikes cutoff idx
    case 2
        %builds and plots timelines for given clusters
colormatrix = colormap('Lines');
load sorted_spikes.mat
spike_num = length(spikes(:,1));
cluster_timeline_2_plot = get(handles.spike_timeline_cluster_to_plot,'String');
cluster_timeline_2_plot = str2num(cluster_timeline_2_plot);


if max(cluster_timeline_2_plot==idx)==1
y=zeros(1,length(ts));
 for i = 1:spike_num
    if idx(i)==cluster_timeline_2_plot

        %plotting
        aaa = find(ts>spike_peak_time(i),1)-1;
        y(aaa:(aaa+length(spikes(i,:))-1))=spikes(i,:);
    end
 end
 plot(ts',y,'Color',colormatrix(cluster_timeline_2_plot,:))
end
if cluster_timeline_2_plot==0
    for cc = 1:max(idx)
    y=zeros(1,length(ts));
    fprintf('Loading cluster %i of %i...\n',cc,max(idx));
    nums2check = find(idx==cc);
    for i = 1:length(nums2check)
    if idx(i)==cc
        %plotting
        try
        aaa = find(ts>spike_peak_time(nums2check(i)),1)-1;
        y(aaa:(aaa+length(spikes(i,:))-1))=spikes(i,:);
        catch lengths_dont_add_up
        end
    end
 end
 hold on
 plot(ts',y,'Color',colormatrix(cc,:))
 clear y
 fprintf('Cluster %i of %i computed...\n',cc,max(idx));
end
end
    case 3
        %plots ISI histograms
        load sorted_spikes
        colormatrix = colormap('Lines');
    for i = 1:max(idx)
    y = spike_peak_time(find(idx==i));
        for ii = 1:length(y)-1
         isi(ii) = y(ii+1)-y(ii);
        end
    subplot(max(idx),1,i)
    hold on
    title(['Cluster ',num2str(i),'   Average Rate: ',num2str(1/mean(isi))])
    isi = sort(isi);
    isi = isi(find(isi>=0));
    isi = isi(find(isi<.3));
    x=hist(isi,100);
    bar(linspace(0,.3,100),x);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',colormatrix(i,:),'EdgeColor',[0,0,0]);
    clear isi
    end
        
    otherwise
end
clear colormatrix
guidata(hObject, handles); %updates the handles



function text_box_for_optional_cluster_override_Callback(hObject, eventdata, handles)
% hObject    handle to text_box_for_optional_cluster_override (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of text_box_for_optional_cluster_override as text
%        str2double(get(hObject,'String')) returns contents of text_box_for_optional_cluster_override as a double
guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function text_box_for_optional_cluster_override_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_box_for_optional_cluster_override (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); %updates the handles
guidata(hObject, handles); %updates the handles

function first_cluster_to_be_combined_Callback(hObject, eventdata, handles)
% hObject    handle to first_cluster_to_be_combined (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of first_cluster_to_be_combined as text
%        str2double(get(hObject,'String')) returns contents of first_cluster_to_be_combined as a double
guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function first_cluster_to_be_combined_CreateFcn(hObject, eventdata, handles)
% hObject    handle to first_cluster_to_be_combined (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); %updates the handles


function second_cluster_to_be_combined_Callback(hObject, eventdata, handles)
% hObject    handle to second_cluster_to_be_combined (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of second_cluster_to_be_combined as text
%        str2double(get(hObject,'String')) returns contents of second_cluster_to_be_combined as a double
guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function second_cluster_to_be_combined_CreateFcn(hObject, eventdata, handles)
% hObject    handle to second_cluster_to_be_combined (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); %updates the handles

% --- Executes on button press in button_to_combine_clusters.
function button_to_combine_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to button_to_combine_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cluster2bcombined = get(handles.first_cluster_to_be_combined,'String');
cluster2bcombined = str2num(cluster2bcombined); %first cluster
if min(cluster2bcombined)>0 
cluster2bcombined = sort(cluster2bcombined); %sorts cluster, uses index of lesser as the new cluster index
combine_clusters(cluster2bcombined);
end
clear cluster2bcombined
set(handles.clusters_to_plot_textbox,'String','0');
set(handles.checkbox_for_plotting_all,'Value',1);
guidata(hObject, handles); %updates the handles
plot_seperated_clusters_button_Callback(hObject, eventdata, handles)

guidata(hObject, handles); %updates the handles



function savefilenametextbox_Callback(hObject, eventdata, handles)
% hObject    handle to savefilenametextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefilenametextbox as text
%        str2double(get(hObject,'String')) returns contents of savefilenametextbox as a double
guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function savefilenametextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefilenametextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); %updates the handles

% --- Executes on button press in bushbutton_to_save_file.
function bushbutton_to_save_file_Callback(hObject, eventdata, handles)
% hObject    handle to bushbutton_to_save_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savefile_name = get(handles.savefilenametextbox,'String');
%savefile_name = strcat(savefile_name,'.mat');
load sorted_spikes.mat
try
save(savefile_name,'idx','spikes','ts','spike_peak_time','interval');
fprintf('Data saved to: ')
fprintf(savefile_name)
fprintf('\n')
clear idx spikes ts savefile_name interval
catch save_file_failed
    fprintf('File could not be saved- missing a variable!\nCheck that "sorted_spikes.mat" contains:\nidx, spikes, ts, spike_peak_time, and interval.\n');
end
guidata(hObject, handles); %updates the handles



function spike_timeline_cluster_to_plot_Callback(hObject, eventdata, handles)
% hObject    handle to spike_timeline_cluster_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spike_timeline_cluster_to_plot as text
%        str2double(get(hObject,'String')) returns contents of spike_timeline_cluster_to_plot as a double
guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function spike_timeline_cluster_to_plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spike_timeline_cluster_to_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); %updates the handles



function cutoff_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to cutoff_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cutoff_textbox as text
%        str2double(get(hObject,'String')) returns contents of cutoff_textbox as a double

guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function cutoff_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cutoff_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

guidata(hObject, handles); %updates the handles


function artifact_catch_textbox_Callback(hObject, eventdata, handles)
% hObject    handle to artifact_catch_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of artifact_catch_textbox as text
%        str2double(get(hObject,'String')) returns contents of artifact_catch_textbox as a double

guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function artifact_catch_textbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to artifact_catch_textbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

guidata(hObject, handles); %updates the handles



function cluster_number_to_be_remove_textd_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_number_to_be_remove_textd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cluster_number_to_be_remove_textd as text
%        str2double(get(hObject,'String')) returns contents of cluster_number_to_be_remove_textd as a double


% --- Executes during object creation, after setting all properties.
function cluster_number_to_be_remove_textd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cluster_number_to_be_remove_textd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cluster_button_to_delete_cluster.
function cluster_button_to_delete_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to cluster_button_to_delete_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cluster2delete = get(handles.cluster_number_to_be_remove_textd,'String');
cluster2delete = str2num(cluster2delete);
remove_cluster(cluster2delete);
set(handles.cluster_number_to_be_remove_textd,'String','0');
set(handles.checkbox_for_plotting_all,'Value',1);
guidata(hObject, handles); %updates the handles
plot_seperated_clusters_button_Callback(hObject, eventdata, handles)

guidata(hObject, handles); %updates the handles


% --- Executes on button press in checkbox_for_plotting_all.
function checkbox_for_plotting_all_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_for_plotting_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_for_plotting_all

guidata(hObject, handles); %updates the handles


% --- Executes on button press in pushbutton_for_batch_clustering.
function pushbutton_for_batch_clustering_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_for_batch_clustering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clc
logf=fopen('logfile.txt','w');fclose(logf); diary logfile.txt;
filestocluster = textread('batchfilelist.txt','%s'); %reads filenames from "batchfilelist.txt" for clustering
for num_of_files = 1:length(filestocluster(:,1))
    try
    set(handles.main_file_name_text,'String',char(filestocluster(num_of_files,:)));
    guidata(hObject, handles); 
    button_main_file_name_Callback(hObject, eventdata, handles)
    if exist('invalidload.mat','file')==0
        load('sorted_spikes.mat','idx')
    if length(idx)>100
    button_to_initiate_clustering_Callback(hObject, eventdata, handles)
    savefilename = strcat('clustered_',filestocluster(num_of_files,:));
    set(handles.savefilenametextbox,'String',char(savefilename));
    guidata(hObject, handles); 
    bushbutton_to_save_file_Callback(hObject, eventdata, handles)
    end
    else
        delete invalidload.mat
        fprintf('Processing next file in list...\n')
        guidata(hObject, handles);
        pause on
        pause(3)
        pause off
    end
    catch failsafeerr
        fprintf('RESET....\n');
    end
    
end
guidata(hObject, handles); %updates the handles



function loadfilenametextbox_Callback(hObject, eventdata, handles)
% hObject    handle to loadfilenametextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadfilenametextbox as text
%        str2double(get(hObject,'String')) returns contents of loadfilenametextbox as a double


% --- Executes during object creation, after setting all properties.
function loadfilenametextbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadfilenametextbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bushbutton_to_load_file.
function bushbutton_to_load_file_Callback(hObject, eventdata, handles)
% hObject    handle to bushbutton_to_load_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadfile_name = get(handles.loadfilenametextbox,'String');
try
load(char(loadfile_name))
save('sorted_spikes.mat','idx','spikes','ts','spike_peak_time','interval');
fprintf('Data Loaded From: ')
fprintf(loadfile_name)
fprintf('\n')
clear idx spikes ts loadfile_name spike_peak_time
catch wrong_load_file_name
    fprintf('File incorrectly formatted!\n');
end
set(handles.clusters_to_plot_textbox,'String','0');
set(handles.checkbox_for_plotting_all,'Value',1);
guidata(hObject, handles); %updates the handles
plot_seperated_clusters_button_Callback(hObject, eventdata, handles)
guidata(hObject, handles); %updates the handles



function spike_length_miliseconds_left_Callback(hObject, eventdata, handles)
% hObject    handle to spike_length_miliseconds_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spike_length_miliseconds_left as text
%        str2double(get(hObject,'String')) returns contents of spike_length_miliseconds_left as a double
guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function spike_length_miliseconds_left_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spike_length_miliseconds_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles); %updates the handles


% --- Executes on button press in button_2_merge_files.
function button_2_merge_files_Callback(hObject, eventdata, handles)
% hObject    handle to button_2_merge_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
filestomerge = textread('mergefilelist.txt','%s'); %reads filenames from "batchfilelist.txt" for clustering
try
for num_of_files = 1:length(filestomerge(:,1))
    load(char(filestomerge(num_of_files,:)))
    if num_of_files ==1
        f_spikes = spikes;
        f_spike_peak_time = spike_peak_time;
        f_idx = idx;
        f_ts = ts;
        next = max(idx);
    else
        if length(f_spikes(1,:))<length(spikes(1,:))
            spikes = spikes(:,1:length(f_spikes(1,:)));
        end
        if length(spikes(1,:))<length(f_spikes(1,:))
            f_spikes = f_spikes(:,1:length(spikes(1,:)));
        end
        f_spikes = [f_spikes;spikes];
        f_spike_peak_time = [f_spike_peak_time;spike_peak_time];
        f_idx = [f_idx;(idx+next)];
        f_ts = [f_ts;ts];
        next = max(f_idx);
    end
end
catch spike_lengths_dissimilar
    fprintf('Spike lengths or intervals dissimilar, stopping\n');
end
spikes = f_spikes; clear f_spikes;
spike_peak_time = f_spike_peak_time; clear f_spike_peak_time;
idx = f_idx; clear f_idx;
ts = f_ts; clear f_ts;
savefile_name = get(handles.merged_final_file,'String');
save(savefile_name,'idx','spikes','ts','spike_peak_time','interval');
fprintf('Files Merged To: ')
fprintf(savefile_name)
fprintf('\n')
set(handles.clusters_to_plot_textbox,'String','0');
set(handles.checkbox_for_plotting_all,'Value',1);
guidata(hObject, handles); %updates the handles
plot_seperated_clusters_button_Callback(hObject, eventdata, handles)

guidata(hObject, handles); %updates the handles



function merged_final_file_Callback(hObject, eventdata, handles)
% hObject    handle to merged_final_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of merged_final_file as text
%        str2double(get(hObject,'String')) returns contents of merged_final_file as a double

guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function merged_final_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to merged_final_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

guidata(hObject, handles); %updates the handles



function old_cluster_name_rename_Callback(hObject, eventdata, handles)
% hObject    handle to old_cluster_name_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of old_cluster_name_rename as text
%        str2double(get(hObject,'String')) returns contents of old_cluster_name_rename as a double


guidata(hObject, handles); %updates the handles
% --- Executes during object creation, after setting all properties.
function old_cluster_name_rename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to old_cluster_name_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

guidata(hObject, handles); %updates the handles


function new_cluster_name_rename_Callback(hObject, eventdata, handles)
% hObject    handle to new_cluster_name_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of new_cluster_name_rename as text
%        str2double(get(hObject,'String')) returns contents of new_cluster_name_rename as a double

guidata(hObject, handles); %updates the handles

% --- Executes during object creation, after setting all properties.
function new_cluster_name_rename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to new_cluster_name_rename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

guidata(hObject, handles); %updates the handles

% --- Executes on button press in pushbutton_to_rename_clusters.
function pushbutton_to_rename_clusters_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_to_rename_clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
old_clus_name = get(handles.old_cluster_name_rename,'String');
old_clus_name = str2num(old_clus_name); %old cluster name
new_clus_name = get(handles.new_cluster_name_rename,'String');
new_clus_name = str2num(new_clus_name); %new
load sorted_spikes.mat
idx(find(idx==old_clus_name))=new_clus_name;
savefile = 'sorted_spikes.mat';
save(savefile,'idx','spikes','ts','spike_peak_time','interval');  
fprintf('Cluster %i renamed to %i\n',old_clus_name,new_clus_name);
set(handles.clusters_to_plot_textbox,'String','0');
set(handles.checkbox_for_plotting_all,'Value',1);
guidata(hObject, handles); %updates the handles
plot_seperated_clusters_button_Callback(hObject, eventdata, handles)
guidata(hObject, handles); %updates the handles


% --- Executes on button press in pushbutton_2_save_as_scr.
function pushbutton_2_save_as_scr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_2_save_as_scr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savefile_name = get(handles.savefilenametextbox,'String');
savefile_name = strcat(savefile_name,'.txt');
load sorted_spikes.mat
fid = fopen(savefile_name,'wt');
interval = ts(2)-ts(1);
%start = find(spikes(1,:)==max(spikes(1,:)));
%spike_start_time = spike_peak_time - interval*start;
ideal_rate = 1/interval;
spike_num = length(idx);
spike_points = length(spikes(1,:));
%header for smr file
fprintf(fid,'%f %i\n',interval,length(spikes(1,:)));
fprintf('Saving File...\n')

for i = 1:length(idx)
fprintf(fid,'%f\t',spike_peak_time(i));
fprintf(fid,'%i\t',idx(i));
fprintf(fid,'%f ',spikes(i,:)); %exports values
fprintf(fid,'\n');
end

fclose(fid);
fprintf('100%% Complete!\n');
disp(['File Saved To: ',savefile_name]);



guidata(hObject, handles); %updates the handles



function spikelengthmiliseconds_right_Callback(hObject, eventdata, handles)
% hObject    handle to spikelengthmiliseconds_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spikelengthmiliseconds_right as text
%        str2double(get(hObject,'String')) returns contents of spikelengthmiliseconds_right as a double


% --- Executes during object creation, after setting all properties.
function spikelengthmiliseconds_right_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spikelengthmiliseconds_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in view_std_aheadoftime.
function view_std_aheadoftime_Callback(hObject, eventdata, handles)
% hObject    handle to view_std_aheadoftime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myfilename = get(handles.main_file_name_text,'String'); %gets filename to be used
cutoff = get(handles.cutoff_textbox,'String'); %gets cutoff
cutoff = str2num(cutoff);
new_data = load(myfilename);
vars = fieldnames(new_data);
load(myfilename,char(vars(1)))
figure
d = char(strcat(vars(1),'.values'));
d = eval(d);
plot(d);
title('Data')
hold on
m = mean(d);
s = std(abs(d));
plot([0,length(d)],[m+cutoff*s,m+cutoff*s],'-r')
plot([0,length(d)],[m-cutoff*s,m-cutoff*s],'-r')
clear filtered d s m amp_thresh spike_num ts spike_peak_time spikes cutoff idx
guidata(hObject, handles); %updates the handles
