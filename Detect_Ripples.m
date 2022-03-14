  
clear; close all;
%%
% Load Filtered file filteredVtraceRp (100-275 Hz) and timestamp or generate below
% Note: LFP should be inverted with respect to
% extracellular convention, i.e., down in the LFP indicates down in spikes
% Input LFP sampling rate
% Input threshold Criterion

ThresholdCriterion = 3; % standard deviations

Fs = 1600;
sampleInterv= 1./Fs;   

[ ts ] = ( 0 : 1 : ( length(filteredVtraceRp) -1 ) ) .* sampleInterv; % generate timestamps  

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect and save ripples from LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% load sleep times
nn = 'Sleeptimes*.mat';
Fnameeeg = dir(nn);
extractnames = extractfield(Fnameeeg,'name');
SWStimes = load(extractnames{:});
SWStimes = SWStimes.SWStimes;   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect: CALCULATE MEAN AND STANDARD DEV. DURING SLEEP LFP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
MarkSleep= zeros(size(SWStimes,1),length(ts));
for s= 1:size(SWStimes,1)
    MarkSleep(s,:)= (ts > SWStimes(s,1) & ts < (SWStimes(s,2)));
end
MarkSleeplogiIdxs = logical(sum(MarkSleep,1)); % LFP
   
sqfiltlfp= filteredVtraceRp.^2;  
m= nanmean(sqfiltlfp(MarkSleeplogiIdxs)); % mean power while animal not moving
sdev= nanstd(sqfiltlfp(MarkSleeplogiIdxs));


% to prevent errors with timing, work with the 'electrodeV' vector with the
% original time relations and NaN when rat not moving
sqfiltlfp= sqfiltlfp.^2;
sqfiltlfp= reshape(sqfiltlfp,1,length(sqfiltlfp)); % make sure it's in the right shape

splog=sqfiltlfp>(m+ (Threshold*sdev)); %%%%% FIRST CRITERION: AT LEAST 3*STD (lfp when quiet) OF FILTERED, SQUARED LFP
pos=find(splog==1);
if isempty(pos)
    ripples{1,1}= 'no ripples in this electrode';
    warndlg('no ripples in this electrode');
else
    difpos=find(diff(pos)>(0.02/sampleInterv)); %%%%%%%% SECOND CRITERION 20ms (at least a difference of 20ms between the end of a ripple and the beginning of the next to be considered separate ripples 
    z=zeros(1,length(difpos)+1); 
    z(1:end-1)=difpos;
    z(end)=length(pos); % now z includes the end position of all ripples
    st= difpos+1; % start positions for all except first ripple
    zst=zeros(1,length(difpos)+1);
    zst(2:end)=st; 
    zst(1)=1;   
    stepripples= [pos(zst);pos(z)]'; % positions in the lfp vector where candidate ripples start and end

%%%%% THIRD CRITERION: eliminate 'ripple' candidates that last less than the sampling interval (one sample)
    duration= stepripples(:,2)-stepripples(:,1);
    stepripples= stepripples((~(duration==0)),:);

   % now determine the start and end based on mean crossing
    below= stepripples(:,1)-(1/sampleInterv); % go back ~1s
    abov= stepripples(:,2)+(1/sampleInterv); % and same above
    tama=size(stepripples);
    tama=tama(1);
    newstepripples=zeros(tama,2);
    for j=1:tama
        if below(j)<0 % if by going back 1s we go below zero
            below(j)=1; % use the first position of the vector instead
        else
        end
        sqfiltlfptemp=sqfiltlfp(below(j):stepripples(j,1));
        for k=length(sqfiltlfptemp):-1:1
           if sqfiltlfptemp(k)< m
               p=abs(length(sqfiltlfptemp)-k);  
               newstepripples(j,1)=stepripples(j,1)-p;
               break
           else
               newstepripples(j,1)=NaN; % if ripple start cannot be found within 1s, don't use it (this occurs when there are NaN in the vector, i.e., animal moving)
           end
        end
    end
         % same for end of ripples
    for j=1:tama
        if abov(j)>length(sqfiltlfp) % if by going above 1s it falls outside the vector
            abov(j)=length(sqfiltlfp); % use the last position of the vector instead
        else
        end
        sqfiltlfptemp=sqfiltlfp(stepripples(j,2):abov(j));
        for k=1:length(sqfiltlfptemp)
            if sqfiltlfptemp(k)< m
                p=k-1;   
                newstepripples(j,2)=stepripples(j,2)+p;
                break
            else
                newstepripples(j,2)=NaN; % if ripple start cannot be found within 1s, don't use it (this occurs when there are NaN in the vector, i.e., animal moving)
            end
        end
    end
    rep= diff(newstepripples(:,1)); 
    newstepripples(rep==0,:)=[]; % eliminate repeated ripples
    incorrect= isnan(newstepripples(:,1)) | isnan(newstepripples(:,2)); % any NaN indicate ripple start/end could not be found because of proximity to a moving period
    newstepripples(incorrect,:)=[];
%%%%%% FOURTH CRITERION:  at least 20ms long  
    duration2= (newstepripples(:,2)-newstepripples(:,1))*sampleInterv; % duration in sec.
    ripplesidx= newstepripples((duration2 > 0.02),:);  % only take those episodes with at least 20ms duration;
    rippletstartIdx= ripplesidx(:,1);
    rippletendIdx= ripplesidx(:,2);
    ripplesIdxs = [rippletstartIdx, rippletendIdx];
end

DurationRipples = (ripplesIdxs(:,2)-ripplesIdxs(:,1))*sampleInterv; % duration in sec. 
rippleTimes = [ts(rippletstartIdx)', ts(rippletendIdx)'];
 
%% SAVE

directorio = pwd;
  
ID1= strfind(directorio,'\');
ID2= ID1(end-1);
name= directorio(ID2:end);
name(regexp(name,'\'))=[];
name= name(find(~isspace(name)));

dirsave = strcat([directorio, '\RippleTimes_' name '_' num2str(Threshold) 'SD']);  
save(dirsave, 'rippleTimes');


