 
clear; close all;
%%
% Load Filtered file (1-4 Hz) and timestamp or generate below
% Note: will detect downward SOs (LFP should be inverted with respect to
% extracellular convention, i.e., down in the LFP indicates down in spikes
% Input LFP sampling rate
% Input threshold Criterion

ThresholdCriterion = 3; % standard deviations

Fs = 1600;
sampleInterv= 1./Fs;   
[ ts ] = ( 0 : 1 : ( length(filteredVtraceSO) -1 ) ) .* sampleInterv; % generate timestamps  

%
 checkfile= ls('FilteredCorticalLFP_1_4_AfterInterpolation_CSC5_*.mat'); %ls(['KcTrough*' epoch '.mat']);  %  ls(['rippletimesCritA*' epoch '.mat']); %ls(['KcTrough*' epoch '.mat']);  %  
 if isempty(checkfile)
     warndlg('File not Found');
 else
     filtLFPSO= load(checkfile); %  ripplestimes= ripplestimes.RPeaks;
     filteredVtraceSO= filtLFPSO.filteredVtraceSO; %ripplestimes.ripplestimes; %
 end      
   
% Check!

figure; plot(ts,filteredVtraceSO,'k');

 
%%
% load sleep times
nn = 'Sleeptimes*.mat';
Fnameeeg = dir(nn);
extractnames = extractfield(Fnameeeg,'name');
SWStimes = load(extractnames{:});
SWStimes = SWStimes.SWStimes;   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect and save SO 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE MEAN AND SDEV DURING SLEEP LFP
  
MarkSleep= zeros(size(SWStimes,1),length(ts));
for s= 1:size(SWStimes,1)
    MarkSleep(s,:)= (ts > SWStimes(s,1) & ts < (SWStimes(s,2)));
end
MarkSleeplogiIdxs = logical(sum(MarkSleep,1)); % LFP
   
relectrodeV= filteredVtraceSO;
relectrodeV(filteredVtraceSO>0)=0; % rectify
SquaredfilteredVtraceSO= relectrodeV.^2;
m= nanmean(SquaredfilteredVtraceSO(MarkSleeplogiIdxs)); % mean power (squared, filtered LFP) during sleep
sdev= nanstd(SquaredfilteredVtraceSO(MarkSleeplogiIdxs));

SquaredfilteredVtraceSO= reshape(SquaredfilteredVtraceSO,1,length(SquaredfilteredVtraceSO));  

splog= SquaredfilteredVtraceSO>(m+(ThresholdCriterion*sdev)); %%%%% FIRST CRITERION: AT LEAST x*STD  OF FILTERED, SQUARED LFP
pos= find(splog==1);

difpos=find(diff(pos)>1);    
z=zeros(1,length(difpos)+1);  
z(1:end-1)=difpos;
z(end)=length(pos); % z includes the end position of all SOs
st= difpos+1; % start positions for all except first 
zst=zeros(1,length(difpos)+1);
zst(2:end)=st; %
zst(1)=1; % zst start position of all SOs
Kon= ts(1,pos(zst)); % get the actual onset time
Kof= ts(1,pos(z));
Kctimes= [Kon;Kof]'; % SOs start and end timestamps
      
  
%% DETECT SO TROUGHS in the filtered LFP

WindowStart= Kctimes(:,1);
WindowEnd= Kctimes(:,2);

TroughSOIDx= zeros(length(WindowStart),1);
TroughKc= zeros(length(WindowStart),1);

goback= 0.02 * 1./sampleInterv;  % go back/forth 20ms
goforth=  0.02 * 1./sampleInterv;  % go back/forth 20ms
TroughAmplitude= zeros(size(WindowStart,1),1);
for klfp= 1:length(WindowStart)
    indi= binsearch(ts,WindowStart(klfp)); % look for the closest time to Kc onset in the time vector  
    indi2= binsearch(ts,WindowEnd(klfp)); % look for the closest time to Kc onset in the time vector  
    if (indi-goback) < 0
        goback= 0;
    elseif (indi2+goforth)> length(filteredVtraceSO)
        goforth= 0;
    end
    lfp= filteredVtraceSO((indi-goback):(indi2+goforth));  % figure; plot(lfp);
    [TT, II]= min(lfp);
    TroughSOIDx(klfp)= indi-goback+(II-1);
    TroughKc(klfp)= ts(indi-goback+(II-1));    
    TroughAmplitude(klfp,1) = abs(TT);
end

 %% SAVE
directorio = pwd;
  
ID1= strfind(directorio,'\');
ID2= ID1(end-1);
name= directorio(ID2:end);
name(regexp(name,'\'))=[];
name= name(find(~isspace(name)));

 
dirsave = strcat([directorio, '\SOTroughtimes_' name '_' num2str(ThresholdCriterion) 'SD']);  
save(dirsave, 'TroughKc');


