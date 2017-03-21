function [spikes spike_peak_time] = obtain_spikes(data,interval,ts,export_the_data,cutoff,milisecleft,milisecright)
% Scott Grimes - Max Planck Cybernetics - 2011
% Obtains spikes from a given channel after filtering and parsing artifacts
signal_begin=round(.5/(milisecleft*1000*interval)); %left side of spike
signal_end=round(.5/(milisecright*1000*interval)); %right side of spike
spike_length=signal_begin+signal_end; %length of spikes
spikes = [zeros(2,spike_length*2)]; 
spike_num=0;
spike_peak_time=0;
t=1;




fprintf('Spike Detection at %3.2fmu...\n',cutoff)
%norm_data = norm_data - m;
%filters using high band pass
% fNorm = (cutoff*s/2);
% [b, a] = butter(10, fNorm, 'high');
% filtered = filtfilt(b, a, norm_data);
% filtered = filtered +m;
% filtered=filtered*pre_norm_max;
clear ATF_check norm_data 
filtered = data;
clear data;

%begin spike detection
sig = std((filtered));
m = mean(filtered);
amp_thresh = cutoff*sig+m;




i=round(signal_begin);
%searches for signals
endcutoff = min([length(ts),length(filtered)]);

while i<endcutoff

    if abs(filtered(i))>amp_thresh
        if i>2*signal_begin && i<length(filtered)-2*signal_end %spike found
        spike_num=spike_num+1; %add 1 to spikes counter
        spikes(spike_num,1:spike_length+1)=[filtered(i-signal_begin:i+signal_end)]; %adds a new row to spikes with the new spike
        spike_peak_time(t,1)=ts(i,1); %gets timestamp of new spike
        t=t+1;
        i=i+signal_end;
       
        end
    end
        i=i+1;
end




if export_the_data==1
savefile = 'sorted_spikes.mat';
save(savefile,'ts','spikes','spike_peak_time','interval')
end
s = sig;
savefile2 = 'ansc1.mat';
save(savefile2,'cutoff','filtered','s','m','amp_thresh','interval')
end