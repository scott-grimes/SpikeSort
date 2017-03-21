function [idx] = combine_clusters(twobcombined)
% Scott Grimes - Max Planck Cybernetics - 2011
% Combines Clusters
load sorted_spikes.mat
for i = 2:length(twobcombined)
 idx(find(idx==twobcombined(i)),1)=twobcombined(1); %set old cluster to new
end
n = histc(idx,1:15);
z = find(n~=0);
for i = 1:length(z)
    idx(find(idx==z(i))) = i;
end
       
       savefile = 'sorted_spikes.mat';
       save(savefile,'idx','spikes','ts','spike_peak_time','interval');
       
fprintf('Clusters ');
fprintf('%i,',twobcombined);
fprintf('\b merged\n');
end