function [LogData] = lin2log_spacing(LinData) 

% The code assumes that your data is in column vectors.

%% Set up parameters
logSpacing = 12;                     % How many sequential points that will have the same spacing before it is doubled.
startAvg = 1;
pointer = 1;
index = 1;
LogData = [];

numColumns = size(LinData,2);
points2avg = startAvg;

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