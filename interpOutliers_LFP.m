function [result] = interpOutliers_LFP(variable,sampleInterv)
% Interpolates small (<500ms) gaps with outliers using linear interpolation
% Marks any outliers (> 10 SD) in the data as
% NaN if longer than 500 miliseconds; if shorter than 500ms: linear interpolation


variable = reshape(variable,1,length(variable)); % make sure vector is horizontal
datacopy.xf = variable;

std_Level = 10; % define level for outliers

% Find outliers >  STD and substitute for NaN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m= mean(datacopy.xf);
s= std(datacopy.xf);

logicdata= (abs(datacopy.xf-m)) > std_Level*s; % find positions of outliers (above or below 6x stand deviation)
len= length(logicdata);
totc=0;
c= 0; % will keep track of number of zeros in a gap
Fs = 1/sampleInterv;
for k= 1:len
    if (logicdata(k)==1) % if there is an outlier
        while logicdata(k)==1
            c= c+1;
            k=k+1;
            if k==len+1
                break
            else
            end
        end
        if c > (Fs/2)  % if more than 500ms of consecutive outlier samples
            %range=[k-c k-1]; 
            datacopy.xf(1,k-c:k-1)= NaN; % do not use the large gaps with outliers
            totc = totc+c; % keep track of the total samples due to outliers
        end
        c=0; % reset the counter
    else
    end 
end
display('Total samples lost because too many consecutive invalid points in the position data:'); % the number of times this message is displayed will = # of gaps too large to be used
totc  % total samples 'lost' due to outliers
  
% Interpolate the small gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logicdata2= find(logicdata==1); % get positions where there are outliers
all= datacopy.xf(logicdata2); % get data in those positions
logicdata3= isnan(all); % determine which ones are large gaps with outliers (they've been substituted for NaN)
unknownx= logicdata2(logicdata3==0); % take the values in the small gaps (these are the x values for which we want to calculate y values in the interpolation)

  
knownx = find(logicdata==0); % known x values -those that are not outliers-
knowny = datacopy.xf(knownx); % corresponding y values

unknowny = interp1(knownx,knowny,unknownx);

datacopy.xf(unknownx)=unknowny;
result= datacopy.xf;


