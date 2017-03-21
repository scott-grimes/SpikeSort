function [idx] = cluster(spikes,optional_cluster_override,pca_yes)
% Scott Grimes - Max Planck Cybernetics - 2011
% Clustering 
L=length(spikes(1,:)); %length of the spike 
spike_num=length(spikes(:,1)); %number of spikes
if pca_yes ==0   %if not using pca...
%wavelet decomp
spikes_wave=0;
for i=1:spike_num
[c,l] = wavedec(spikes(i,:),5,'db5');
if i==1
    spikes_wave=c;
else
spikes_wave=[spikes_wave;c];
end
end
fprintf('Wavelet Decomp Complete...\n')
clear c
%ends wavelet decomp
idx = zeros(spike_num,1);
%normalizing spikes
for i = 1:spike_num
    spikes(i,1:L) = spikes(i,:)/norm(spikes(i,:));
end

 
if optional_cluster_override==0 %auto detect cluster number
idx = kmeans([spikes_wave],15,'EmptyAction','drop');
 
k = cluster_consol([spikes_wave],idx);
idx = kmeans([spikes_wave],k,'EmptyAction','drop');
fprintf('Using %i-Clusters\n',k);
else %manual cluster number
    idx = kmeans([spikes_wave],optional_cluster_override,'EmptyAction','drop');
end




else %using PCA 
[coef,score,latent,tsquare]=princomp(spikes);
fprintf('PC loadings calculated...\n'); 
comp = cumsum(latent)./sum(latent); 
k = find(comp>=.98,1);
if optional_cluster_override==0
idx = kmeans(score(:,1:k),15,'EmptyAction','drop');
k = cluster_consol(score(:,1:k),idx);
idx = kmeans(score(:,1:k),k,'EmptyAction','drop');
fprintf('Using %i-Clusters\n',k);
clear coef
else
idx = kmeans(score(:,1:k),optional_cluster_override,'EmptyAction','drop'); %manual cluster number   
end
    








end

fprintf('Sorting Completed - %i Clusters\n',max(idx))

