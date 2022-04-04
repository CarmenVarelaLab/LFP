function [downsampledLFP] = downsampleLFP(RawLFP,sr,dfactor)
% Downsample the LFP by a given (integer) factor using 'resample'...
...'resample' applies an FIR Antialiasing Lowpass Filter...
    ...to x and compensates for the delay introduced by the filter...
    ...Also saves the downsampledLFP to the local directory
    
if isempty(dfactor)
    dfactor = 20;
    warndlg('No downsampling factor provided. Using default= 20');
else
end

RawLFP = double(RawLFP);
RawLFP = reshape(RawLFP,1,length(RawLFP)); % make sure vector is horizontal
tot = length(RawLFP);

% To avoid running out of memory issues will work with chunks:
% first reshape the LFP vector into a matrix
if rem(tot, 2) == 0 % make even if not
else
    LastSample = RawLFP(end);
    RawLFP = RawLFP(1:end-1);  % drop last sample
end

if tot > 20000
    K = 1:tot;
    D = K(rem(tot,K)==0); % find divisors of tot
    D2 = D(rem(D,dfactor) ==0); % elements in D that will give an integer value when downsampled by dfactor
    
    [fd,~] = findNearest(D2,10000); % find a value near 10000 that can be used to reshape the LFP vector into a matrix
    ReshapedLFP = reshape(RawLFP,fd,[]);
    
    %%  Downsample
    
    downsampled = zeros(size(ReshapedLFP,1)./dfactor,size(ReshapedLFP,2));
    for d = 1:size(downsampled,2)
        temp = ReshapedLFP(:,d);
        downsampled(:,d)= resample(temp',1,dfactor);
    end
    downsampledLFP = reshape(downsampled,1,[]);
    
    %
    % sampleInterv = 1 ./ 25000;
    % tot = length( RawLFP );
    % [ ts ] = ( 0 : 1 : ( tot -1 ) ) .* sampleInterv; % generate timestamps with old sampling rate
    %
    % ts(end)
    % ts2 = linspace(ts(1),ts(end),length(downsampledLFP));
    %
    % figure; plot(ts,RawLFP,'k');
    % hold on; plot(ts2,downsampledLFP,'r');
    if rem(tot, 2) == 0 % make even if not
    else
        downsampledLFP = [downsampledLFP,LastSample]; % keep the same number of samples as the original       
    end
     
    %% save to current directory
    directorio = pwd;
    
    ID1= strfind(directorio,'\');
    ID2= ID1(end-1);
    name= directorio(ID2:end);
    name(regexp(name,'\'))=[];
    name= name(find(~isspace(name)));
    
    NewSR = sr./dfactor;
    
    dirsave = strcat([directorio, '\LFP_Downsampled' name '_' num2str(NewSR) 'Hz']);
    save(dirsave, 'downsampledLFP');
    
else
    warndlg('Does not run for LFPs shorter than 20000 samples');
end

