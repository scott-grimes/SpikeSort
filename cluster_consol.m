function y = cluster_consol(spikes,idx)
% Scott Grimes - Max Planck Cybernetics - 2011
% consolidates Clusters
L = length(spikes(1,:)); %length of spikes
idx_old=idx; %saves original spike index
k = max(idx); %number of clusters

varr=ones(k,k); %var matrix

for i = 1:k
    c(i,1:L) = group_clusters(spikes,idx,i); %groups the spikes by their cluster index
end
fprintf('Clusters Grouped, generating sim matrix...\n');
    k = max(idx);
    varr=zeros(k,k);
for center = 1:k  %gets varriances
    for test = 1:k
        if test>center
            
        varr(center,test)= dot(c(center,1:L),c(test,1:L))/(sum(c(center,1:L).^2)^.5*sum(c(test,1:L).^2)^.5); %sets var_nm to the var between clusters n, m
        else
            varr(center,test)=0; %lower half of matrix left blank
        end
    end
end
kval = max(idx_old); %initial k guess > final k
maxvarval = max(max(varr)); %selects clusters with greatest similarity

while length(varr(1,:))>2 %combines clusters with greatest similarity until only 2 clusters remain
    
       [i throwaway]=  find(max(max(varr))==varr,1,'first'); %gets cell with max var
       idx(find(idx==throwaway),1)=i;
       for zz = 1:length(idx)
           if idx(zz)>throwaway
               idx(zz)=idx(zz)-1;
           end
       end
       k = max(idx); %new max k
       fprintf('Generating new similarities: %i to %i...\n',k+1,k);
       
       %calc'ing new var's for new row...
       c = c([1:throwaway-1,throwaway+1:k+1],1:L);
       c(i,1:L) = group_clusters(spikes,idx,i)';
       varr = varr([1:throwaway-1,throwaway+1:k+1],[1:throwaway-1,throwaway+1:k+1]);
       
       for newvals = 1:i-1
           varr(newvals,i) = dot(c(newvals,1:L),c(i,1:L))/(sum(c(newvals,1:L).^2)^.5*sum(c(i,1:L).^2)^.5);
       end
       
       for newvals = i+1:k
           varr(i,newvals) = dot(c(newvals,1:L),c(i,1:L))/(sum(c(newvals,1:L).^2)^.5*sum(c(i,1:L).^2)^.5);
       end
       
       kval = [kval k];
       maxvarval = [maxvarval max(max(varr))];
end


ydd = diff(maxvarval,2); %second prime of var's that were combined
y = 1;
for i = 1:length(ydd)-1;
if i==1
ycheck = ydd(i);
else
ycheck = [ycheck ydd(i)];
m = mean(ycheck);
s = std(ycheck);
if abs(ydd(i+1))>m+3*s %checks for outliers
y = i+1;
break
end
end
end
k = kval(kval(1)-y); %determines k final value
if k<2
    k=2; %if k<2 set k=2 (cannot have one cluster)
end
%figure
%plot(ydd,'*')
y = k;
end