% Time = [];
% Ydata = [];
% positive_time = [];
% positive_time_data = [];
% T0_Time = [];
% TimeAbove0 = [];
% DataAbove0 = [];
% LogTime = [];
% intensities = [];
% intensityatT0 = [];
% FitX = [];
% FitY = [];
% FitFunctions = {};
% t50_time_evaluator = [];

% Assign variables
Time = LogTimeArray(:,1);
Ydata = LogAbsArray(:,3);

% Perfom filtering of data
positive_time = Time(Time > 0);
positive_time_data = Ydata(Time > 0, :);

DataTime = positive_time;
DataTable = positive_time_data;

T0_Time = DataTime(DataTable == max(DataTable));

% Load the time
TimeAbove0 = DataTime(DataTime >= T0_Time);
DataAbove0 = DataTable(DataTime >= T0_Time);

% Defining the fit functionfor eval
linfittype = fittype('a/(1+b*10^x)^c + d');

% Define the fit options for the line
fopt = fitoptions('Method','NonLinearLeastSquares');

% These are the values correspondiong to a,b,c and d
fopt.Lower = [-inf -inf 0 -inf];
fopt.Upper = [inf inf 1 inf];

% Obtain the logarithm of all time and y above 0
LogLogTime = log10(TimeAbove0);
data = DataAbove0;

% Prepare data of same size (same columns)
[rowsize,colsize]=size(data);
% It's less memory exhaustive to use loop
% We know data size to expect
FitX = zeros(rowsize,1);
% Keeps track of all fittedY
FitY = zeros(rowsize, colsize);
FitFunctions = cell(colsize,1);

for column=1:1:colsize
    % prepare fit data for current kinetic data
    [FitX, FitY(:,column)] = prepareCurveData(LogLogTime, data(:,column));
    % set starting intensity to 0 to be able to handle both positive and negative initial amplitudes
    % set infinite baseline to the last index of the fitY data
    fopt.StartPoint = [0 4 0.5 FitY(end, column)];
    % Define range of data which to exclude (now redundant - data only has values for time > t0)
    % Extremeties don't work with this line enabled
    %fopt.Exclude = app.FitX < StartFitTime;
    FitFunctions{column} = fit(FitX, FitY(:,column), linfittype, fopt);
end

% Get T0 intensity from fit function
intensityatT0 = feval(FitFunctions{1,1}, log10(T0_Time));
% Get rest of intensities as fit
intensities = feval(FitFunctions{1,1},FitX);

% Obtain t50% intensity
t50_intensity = intensityatT0/2;
% Define function for obtaining the time using given intensity
t50_time_evaluator = @(x)t50_intensity - FitFunctions{1,1}(x);
% Obtain the fitted time from fit function
logt50_time = fzero(t50_time_evaluator, log10(T0_Time));
% Convert time to base 10
t50_time = 10^logt50_time;
disp(t50_time)

f = figure;
ax = axes(f);

plot(ax,DataTime,DataTable);
line(ax,TimeAbove0,intensities,'LineStyle','-','Color','r')
hold on
xline(ax,t50_time);
hold off
ax.XScale = 'log';