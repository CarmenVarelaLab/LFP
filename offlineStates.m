function SleepPeriods = offlineStates(delta,ripple,window,step,sr)
% INPUTs 
%  delta    vector with delta filtered LFP
%  ripple     vector with ripple filtered LFP
%  window    length (in secs) of the window to average the LFP
%  step     length (secs) of the steps for the moving average
%  sr       sampling rate (Hz)
%
% OUTPUT 
% 	SleepPeriods = start and end times of putative quiet wakefulness/sleep
% 	periods

deltasq= delta.^2;   
ripplesq= ripple.^2;


win= (window*sr)-1;  
stepsize= step.* sr; % moving in steps of 1s  (600samples)

steps= 1:stepsize:(length(deltasq)-1); % steps of 1s

% corresponding timestamps for the new vectors
Avtimest= steps./sr;

AvDelta= zeros(1,length(steps));
AvRipple= zeros(1,length(steps));
c=1;
for i= steps
    if i+win>(length(deltasq))
        AvDelta(c)= NaN;
        AvRipple(c)= NaN;
        break
    else
        AvDelta(c)= mean(deltasq(i:i+win));
        AvRipple(c)= mean(ripplesq(i:i+win));
    end
    c=c+1;
end


mD= nanmean(AvDelta);
mR= nanmean(AvRipple);

Power= [AvDelta; AvRipple];

SWS= (Power(1,:)> mD & Power(2,:)> mR);  %%%%%%% FIRST CRITERION: BOTH RIPPLE AND DELTA POWER ABOVE AVERAGE

pos=find(SWS==1);

difpos=find(diff(pos)>(120 .* step)); %%%%%%%% SECOND CRITERION: at least a difference of 120sec between SWS periods to be considered separate (both average delta and ripple power are now sampled in 1sec steps

z=zeros(1,length(difpos)+1);  
z(1:end-1)=difpos;
z(end)=length(pos); %  z includes the end position of all putative sleep periods
st= difpos+1; % start positions for all except first sleep period

zst=zeros(1,length(difpos)+1);
zst(2:end)=st; %
zst(1)=1; % zst start position for first
stepSWS= [pos(zst);pos(z)]'; % positions in the original vector where candidate SWS periods start and end

SWSstart= Avtimest(1,stepSWS(:,1));
SWStend= Avtimest(1,stepSWS(:,2));
SWStimes= [SWSstart; SWStend]';


dur=SWStimes(:,2)-SWStimes(:,1);

SleepPeriods= SWStimes(dur>120,:);

