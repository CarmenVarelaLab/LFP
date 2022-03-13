

% Note: interpolates LFP outlier points before filtering: points  > 6 standard deviations; if outliers for 500ms or more--- they are substituted by NaNs, and the total number of samples is displayed
%
clear; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input file information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = {'Enter Cx CSC channel','Enter HC CSC channel','Enter Pos File Name','Enter Events File Name'};
dlg_title = 'Enter channels for analysis';
num_lines = 1;
def = {'CSC5_LFP','CSC9_LFP','VT1','Events'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
Cx= answer{1,1};
HC= answer{2,1};
Pos= answer{3,1};
Evnt = answer{4,1};

CxLFPFile = strcat([Cx '.ncs']); % Cx LFP name
HCLFPFile = strcat([HC '.ncs']); % HC LFP name
PosFile = strcat([Pos '.nvt']); % Pos File name
EventsFile = strcat([Evnt '.nev']); % Events File name

% Conversion factors for velocity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prompt = {'Enter cms X side','Enter cms Y side','Enter pixel range X side','Enter pixel range Y side','Camera s.r.'};
dlg_title = 'Enter cms/Pixel correspondence for Vel';
num_lines = 1;
def = {'35','38','240','262','30'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
cmsSideX= str2double(answer{1,1});
cmsSideY= str2double(answer{2,1});
PixelsSideX= str2double(answer{3,1});
PixelsSideY = str2double(answer{4,1});
srCamera = str2double(answer{5,1});

sampleintCam = 1/srCamera; % camera sampling interval

cmPpixelX= cmsSideX/PixelsSideX; % range of X pixels with camera lens at a particular height from sleep box floor 
cmPpixelY= cmsSideY/PixelsSideY; % range of Y pixels 
cmfac = mean([cmPpixelX, cmPpixelY]);

% parameters for sleep/quiet classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Enter velocity Threshold to define quiet wake/sleep','Enter number of samples for smoothing'};
dlg_title = 'Enter cms/Pixel correspondence for Vel';
num_lines = 1;
def = {'2','15'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

Threshold_Vel= str2double(answer{1,1}); % below 2cms/s --- quiet/sleep
MovMeanSamples= str2double(answer{2,1}); % number of samples to smooth rat coordinates for 0.5 second moving mean


% genearte name to save files later
directorio = pwd;

ID1= strfind(directorio,'\');
ID2= ID1(end-1);
name= directorio(ID2:end);
name(regexp(name,'\'))=[];
name= name(find(~isspace(name)));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract LFP and timestamp data from Neuralynx CSC files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Cortical LFP
[ timestamps, sampleFrequencies, CxLFP, header ]= Nlx2MatCSC(...
    CxLFPFile,...
    [ 1 0 1 0 1 ],...
    1,...
    1,...
    [ ] );
%% to load just part of the file
% 
% [ timestamps, sampleFrequencies, CxLFP, header ]= Nlx2MatCSC(...
%     CxLFPFile,...
%     [ 1 0 1 0 1 ],...
%     1,...
%     4,...
%     [1632912364203345 (1632912364203345+200000) ] );


%%
  
is_present = not( cellfun( @isempty, strfind( header, 'ADBitVolts' ) ) );

ADBitVolt = header{is_present ==1};
ADBitVolt(1:12) = [];
ADBitVolt = str2num(ADBitVolt);

CxLFP = CxLFP(:)' .* ADBitVolt .* 1000; %  miliVolts
 
if isempty(HC)  % if working only with Cx skip
else    
    % HC LFP
    [ timestamps, sampleFrequencies, HCLFP, header ]= Nlx2MatCSC(...
        HCLFPFile,...
        [ 1 0 1 0 1 ],...
        1,...
        1,...
        [ ] );
    
    HCLFP = HCLFP(:)'.* ADBitVolt .* 1000; % in miliVolts
end

% generate timestamps
sampleInterv = 1 ./ sampleFrequencies ( 1 );
tot = length( CxLFP );
[ ts ] = ( 0 : 1 : ( tot -1 ) ) .* sampleInterv; % generate timestamps


%% interpolate outliers > 6 standard deviations; if outliers for 500ms or more--- they are substituted by NaNs, and the total number of samples is displayed

oldCxLFP = CxLFP;
oldHCLFP = HCLFP;

CxLFP = interpOutliers_LFP(CxLFP,sampleInterv);
HCLFP = interpOutliers_LFP(HCLFP,sampleInterv);

% Check LFPs look ok
figure; plot(ts, oldHCLFP,'k');
hold on; plot(ts,HCLFP,'b');

figure; plot(ts, oldCxLFP,'k');
hold on; plot(ts,CxLFP,'b');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter Cortical LFP in SO and spindle bands  & HC LFP in Ripple Band
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

Fs= 1/sampleInterv; % sampling rate Hz

% SO band
%%%%%%%%%%
lfreqSO= 1;  % lower limit
hfreqSO= 4;  % upper limit

band=[lfreqSO hfreqSO];
bandpass= hfreqSO-lfreqSO;

n= ceil(6*(Fs/bandpass));  % minimum filter order when using a blackman window = 6 * (Fs/main-lobe-bandwidth) (from 'filter_design' pdf, page 5, references Mitra p.535
%   the n in the Matlab function fir1() is the filter order, i.e. you get a
%   vector with n+1 elements as a result (so n+1 is your filter length =
%   number of taps) http://dsp.stackexchange.com/questions/8685/filter-order-vs-number-of-taps-vs-number-of-coefficients
% make sure n is even
if rem(n,2)==0
    n=n;
else % if not...
    n=n+1;
end

b = fir1(n, 2*band./Fs,blackman(n+1));

filteredVtraceSO = filtfilt(b,1,CxLFP);


% spindle band
%%%%%%%%%%%%%%%%%%

lfreqSp= 6;  % lower limit
hfreqSp= 14;  % upper limit

band=[lfreqSp hfreqSp];
bandpass= hfreqSp-lfreqSp;

n= ceil(6*(Fs/bandpass));
if rem(n,2)==0
    n=n;
else % if not...
    n=n+1;
end

b = fir1(n, 2*band./Fs,blackman(n+1));

filteredVtraceSP = filtfilt(b,1,CxLFP);

if isempty(HC)
else
    % ripple band
    %%%%%%%%%%%%%%%%

    lfreqRp= 100;  % lower limit
    hfreqRp= 275;  % upper limit

    band=[lfreqRp hfreqRp];
    bandpass= hfreqRp-lfreqRp;

    n= ceil(6*(Fs/bandpass));
    if rem(n,2)==0
        n=n;
    else % if not...
        n=n+1;
    end

    b = fir1(n, 2*band./Fs,blackman(n+1));

    filteredVtraceRp = filtfilt(b,1,HCLFP);
end
% 

%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save Raw and Filtered LFPs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd ..
mkdir Analysis;
cd('Analysis');

directorio = pwd;

dirsave = strcat([directorio, '\RawCorticalLFP_AfterInterpolation' Cx '_' name]); 
save(dirsave, 'CxLFP');

dirsave = strcat([directorio, '\Timestamp_' name]); 
save(dirsave, 'ts');

% filtered signals
dirsave = strcat([directorio, '\FilteredCorticalLFP_1_4_AfterInterpolation_' Cx '_' name]); 
save(dirsave, 'filteredVtraceSO');

dirsave = strcat([directorio, '\FilteredCorticalLFP_6_14_AfterInterpolation_' Cx '_' name]); 
save(dirsave, 'filteredVtraceSP');

if isempty(HC)
else    
    dirsave = strcat([directorio, '\RawHippocampalLFP_AfterInterpolation_' HC '_' name]);
    save(dirsave, 'HCLFP');
    
    dirsave = strcat([directorio, '\FilteredHippocLFP_100_275_AfterInterpolation_' HC '_' name]);
    save(dirsave, 'filteredVtraceRp');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process position info, Calculate Velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[TimeStampsPos, X, Y] = Nlx2MatVT(PosFile, [1 1 1 0 0 0],0,1,1);
totPos = length(X);

tsPos = (0 : 1 : (totPos-1)).*sampleintCam;

% interpolate zeros & outliers (above or below 6 standard deviations) that
% occur for less than 9samples (9*0.033= 0.297secs); if longer, it marks
% the gaps with NaNs and keeps track of total time with NaNs

Xinterp1 = interpinval_Coordinates(X);
Yinterp1 = interpinval_Coordinates(Y);


Xinterp = movmean(Xinterp1,MovMeanSamples,'omitnan');
Yinterp = movmean(Yinterp1,MovMeanSamples,'omitnan');

XYinterp = [Xinterp;Yinterp];
 
figure; plot(tsPos,Y,'.k');
hold on; plot(tsPos,Yinterp,'.r');
figure;  plot(tsPos,Xinterp,'.b');

% Velocity

Vel= gradient(Xinterp, sampleintCam) + i*gradient(Yinterp, sampleintCam); % instantaneous velocity (derivative with dt=33msec
fvel= abs(Vel); % Instantaneous speed is the magnitude of instantaneous velocity.
velcms= fvel.* cmfac; %  vel magnitude in cm/s

%  figure; plot(velcms,'r');

Velocity.VelCmsPSec = velcms;
Velocity.Timestamp = tsPos;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dirsave = strcat([directorio, '\XYcoordinates_' name]); 
save(dirsave, 'XYinterp');

dirsave = strcat([directorio, '\VelocityRaw_' name]); 
save(dirsave, 'Velocity');
