function y = group_clusters(spikes1,idx1,k)
    tempidx = ((idx1==k)*k==0);
    for i=1:length(idx1)
    temp(i,:) = spikes1(i,:)*tempidx(i);
    end
    y = sum(temp)/sum(tempidx);
end