function [idx] = remove(getridof)
% Scott Grimes - Max Planck Cybernetics - 2011
% Removes Clusters
load sorted_spikes.mat

for i = 1:length(getridof)
idx(find(idx==getridof(i)))=0;  %sets cluster to be removed indexes to zero
spikes = spikes(find(idx~=0),:); %removes spikes with index of 0 
idx = idx(find(idx~=0),:); % removes indexes with value of 0
spike_peak_time = spike_peak_time(find(idx~=0),:); %removes timestamps with index value of 0
end
n = histc(idx,1:15);
z = find(n~=0);
for i = 1:length(z)
    idx(find(idx==z(i))) = i;
end
       
       savefile = 'sorted_spikes.mat';
       save(savefile,'idx','spikes','ts','spike_peak_time','interval');
end