function [FR,SWB]=FR_SWB(spiketime,TT)
dt=0.02;
% firing rate
     if spiketime==0
        FR=0;
     else
     FR=length(spiketime)./(TT/10^3); % mean frequency
     end

% spiketime in samples (ms/dt)
% number of spikes

         if isempty(spiketime)==1
        spikesN=0;
         else
    spikesN=length(spiketime);
         end

% for %SWB
% burst detection

            if spiketime==0
            brst=[];
            else
            brst=buda_detect_bursts_canonical(spiketime,80/dt,160/dt,2); % classical 80/160 rule            
            end
% ISI
            isi=diff(spiketime); 

% SWB
        if isempty(brst)==1
        totalspikesbrst=0;
        else
        totalspikesbrst=sum(brst.nSp); % total n of spikes in burst
        end
        totalspikes=length(isi)+1; % total n of spikes
        SWB=totalspikesbrst./totalspikes; % standart burst measure