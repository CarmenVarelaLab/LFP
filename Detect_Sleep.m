
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLEEP DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD: LFP timestamp (ts);
% delta-filtered LFP (filteredVtraceSO);
% ripple-filtered LFP (filteredVtraceRp);
% velocity
% needs toPPT to transfer figures to powerpoint: https://www.mathworks.com/matlabcentral/fileexchange/44851-jrichter24-toppt
% Note:  saves sleep BOUTS THAT ARE AT LEAST 1 MIN LONG
 

Fs = 1600;

deltabandP= filteredVtraceSO.^2;
deltabandH= abs(hilbert(deltabandP));
   
  
if isempty(HC)
else
    ripplbandP= filteredVtraceRp.^2;
    ripplbandH= abs(hilbert(ripplbandP));
end
  
% calculate average with moving window

win= (2*Fs)-1; % 2sec window
stepsize= Fs; % moving in steps of 1s   
steps= 1:stepsize:(length(deltabandH)-1); % steps of 1sec

AvDelta= zeros(1,length(steps));
AvRipple= zeros(1,length(steps));
c=1;
for i= steps
    if i+win>(length(deltabandH))
        break
    else
        AvDelta(c)= mean(deltabandH(i:i+win));
        if isempty(HC)
        else          
            AvRipple(c)= mean(ripplbandH(i:i+win));              
        end
    end
    c=c+1;
end

% corresponding timestamps for the new vectors

Avtimest= steps./Fs;

% assign SWS if delta and ripple > mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mD= mean(AvDelta);
sdDelta = std(AvDelta);

if isempty(HC)
    SWS= AvDelta > (mD);   %%%%%%% FIRST CRITERION: DELTA POWER ABOVE AVERAGE
else    
    mR= mean(AvRipple);    
    Power= [AvDelta; AvRipple];     %%%%%%% IF HC LFP AVAILABLE-- FIRST CRITERION: DELTA and RIPPLE POWER ABOVE AVERAGE
    SWS= (Power(1,:)> mD & Power(2,:)> mR);  %%%%%%% SLEEP CRITERION: BOTH RIPPLE AND DELTA POWER ABOVE AVERAGE
end

%%%%%%%%
pos=find(SWS==1);
  
difpos=find(diff(pos)>(10)); 
z=zeros(1,length(difpos)+1); 
z(1:end-1)=difpos;
z(end)=length(pos); %  z includes the end position of all candidate sleep bouts
st= difpos+1; % start positions for all except first bout
zst=zeros(1,length(difpos)+1);
zst(2:end)=st; %
zst(1)=1;  
stepSWS= [pos(zst);pos(z)]'; % positions in the original vector where candidate SWS periods start and end

incorrect= isnan(stepSWS(:,1)) | isnan(stepSWS(:,2));   
stepSWS(incorrect,:)=[];

SWSstart= Avtimest(1,stepSWS(:,1));
SWStend= Avtimest(1,stepSWS(:,2));
SWStimes= [SWSstart; SWStend]';


dur=SWStimes(:,2)-SWStimes(:,1);
SWStimes= SWStimes(dur>60,:);  %%%%% KEEP BOUTS THAT ARE AT LEAST 1 MIN LONG

%%
% plots
scrsz=get(0,'screensize');
f1 = figure('Position',[80 80 scrsz(3) scrsz(4)/1.2]);
hold on;
timest= Avtimest;
% if both cortex and HC LFPs were used
if isempty(HC)    
    % plot Av delta power
    subplot(2,2,[1 2]);
    plot(timest,AvDelta,'k');
    hold on;
    line([0 timest(end)],[mD mD], 'color','k','linestyle','--','linewidth',2);
    set(gcf, 'color', 'white');
    set(gca, 'Box','off','fontsize',12, 'fontweight', 'bold');
    %set(gca, 'xcolor','w'); % no xaxis
    set(gca,'linewidth',1.2);
    set(gca,'TickDir','out');
    ylabel('Av Delta Power');
    xlim([0 timest(end)]);
    
    subplot(2,2,[3 4]);
    plot(tsPos,velcms,'b');
    hold off;
    set(gcf, 'color', 'white');
    set(gca, 'Box','off','fontsize',12, 'fontweight', 'bold');
    %set(gca, 'xcolor','w'); % no xaxis
    set(gca,'linewidth',2.1);
    set(gca,'TickDir','out');
    ylabel('Vel (cm/s)');
    xlim([0 tsPos(end)]);
    
    ylim([0 50]);
else        
    subplot(3,2,[1 2]);
    plot(timest,AvRipple,'k');
    hold on;
    line([timest(1) timest(end)],[mR mR], 'color','k','linestyle','--','linewidth',2);
    set(gcf, 'color', 'white');
    set(gca, 'Box','off','fontsize',12, 'fontweight', 'bold');
    set(gca,'linewidth',1.2);
    set(gca,'TickDir','out');
    ylabel('Av Rpl Power');
    xlim([0 timest(end)]);    
    
    % add SWS periods in same subplot
    hold on;
    
    SWSplot= SWStimes;
    ta=size(SWSplot);
    ta=ta(1);
    for i=1:ta
        line([SWSplot(i,1) SWSplot(i,2)], [0.0007 0.0007],'color',[0.5843 0.8157 0.9882],'linestyle','-','linewidth',4);
    end
    ylim([0 0.0008]);
    namet =name;  
    namet(regexp(namet,'_'))= ' ';
    title(namet,'VerticalAlignment','bottom','HorizontalAlignment','left');
    
    subplot(3,2,[3 4]);
    plot(timest,AvDelta,'k');
    hold on;
    line([0 timest(end)],[mD mD], 'color','k','linestyle','--','linewidth',2);
        
    hold on;
    
    SWSplot= SWStimes;
    ta=size(SWSplot);
    ta=ta(1);
    for i=1:ta
        line([SWSplot(i,1) SWSplot(i,2)], [0.0035 0.0035],'color',[0.5843 0.8157 0.9882],'linestyle','-','linewidth',4);
    end
    
    set(gcf, 'color', 'white');
    set(gca, 'Box','off','fontsize',12, 'fontweight', 'bold');
    %set(gca, 'xcolor','w'); % no xaxis
    set(gca,'linewidth',1.2);
    set(gca,'TickDir','out');
    ylabel('Av Delta Power');
    xlim([0 timest(end)]);
    
    subplot(3,2,[5 6]);    
    plot(tsPos,velcms,'b');
    hold off;
    set(gcf, 'color', 'white');
    set(gca, 'Box','off','fontsize',15, 'fontweight', 'bold');
    set(gca,'linewidth',2.1);
    set(gca,'TickDir','out');
    ylabel('Vel (cm/s)');
    xlim([0 tsPos(end)]);    
   % ylim([0 50]);    
end

% Send figure to powerpoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toPPT(f1,'Width%',80,'SlideNumber','append','exportFormatType','bmp');

% save SWS times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% genearte name to save files later
directorio = pwd;

ID1= strfind(directorio,'\');
ID2= ID1(end-1);
name= directorio(ID2:end);
name(regexp(name,'\'))=[];
name= name(find(~isspace(name)));

dirsave= strcat([directorio '\Sleeptimes' name]);
save (dirsave, 'SWStimes');

%%  BASIC METRICS
% Total sleep
TotSleepSecs = sum(SWStimes(:,2)-SWStimes(:,1));
% Fraction of time sleeping
FracSleep = (TotSleepSecs./ts(end)).*100;




