function [LogTime, LogData] = SpectraBuilder_lin2log_TAS(LinTime,LinData)

% The code assumes that your data is in column vectors.

%% Set up parameters
logSpacing = 8;                     % How many sequential points that will have the same spacing before it is doubled.
startAvg = 1;
LogData = [];

numColumns = size(LinData,2);

%% Find time = 0 index

ZeroIndex = find(abs(LinTime) == min(abs(LinTime)));

%% Calculate for the time data.

%Calculate for time < 0.
points2avg = startAvg;
pointer = ZeroIndex;
index = 1;

while pointer-(logSpacing*points2avg) > 1
    
    for repeat = 1:1:logSpacing
        
        LogTime(index) = mean(LinTime(pointer-points2avg:pointer));
                
        pointer = pointer - points2avg;
        index = index +1;
    end
    
    points2avg = 2*points2avg;
end

LogTime = fliplr(LogTime); % To align time from most negative up to zero.

%Calculate for time time > 0.

points2avg = startAvg;
pointer = ZeroIndex;

while pointer+(logSpacing*points2avg) < length(LinTime)
    
    for repeat = 1:1:logSpacing
        
        LogTime(index) = mean(LinTime(pointer:pointer+points2avg));
                
        pointer = pointer + points2avg;
        index = index +1;
    end
    
    points2avg = 2*points2avg;
end

LogTime = LogTime';

%% Calculate for the y data

%Calculate for time < 0.
points2avg = startAvg;
pointer = ZeroIndex;
index = 1;

while pointer-(logSpacing*points2avg) > 1
    
    for repeat = 1:1:logSpacing
        
        for column = 1:1:numColumns
        LogData(index,column) = mean(LinData(pointer-points2avg:pointer,column));
        end
        
        pointer = pointer - points2avg;
        index = index +1;
    end
    
    points2avg = 2*points2avg;
end

LogData = flipud(LogData); % To align for time from most negative up to zero.

%Calculate for time < 0.
points2avg = startAvg;
pointer = ZeroIndex;

while pointer+(logSpacing*points2avg) < length(LinData)
    
    for repeat = 1:1:logSpacing
        
        for column = 1:1:numColumns
            LogData(index,column) = mean(LinData(pointer:pointer+points2avg,column));
        end
        
        pointer = pointer + points2avg;
        index = index +1;
    end
    
    points2avg = 2*points2avg;
end