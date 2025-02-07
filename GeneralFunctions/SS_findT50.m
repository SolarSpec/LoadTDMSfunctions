function t50 = SS_findT50(cfitObject,t0,ncell,BaselineIndex)
% Find the t50% from the input time 0, using the cfitObject. Look for the
% ncell of the cfit cell array outputted by the global fitting app.
% BaselineIndex is the index of which coefficient is the baseline offset.

cFitFunction = cfitObject{ncell}

FitCoeffs = coeffvalues(cFitFunction);
BaselineCoeff = FitCoeffs(BaselineIndex); 

% Get T0 intensity from fit function
intensityatT0 = feval(cFitFunction, t0);

% Obtain t50% intensity of decay (i.e., consider the baseline offset)
t50_intensity = (intensityatT0+BaselineCoeff)/2;

% Define function for obtaining the time using the given intensity
t50_time_evaluator = @(x)t50_intensity - cFitFunction(x);

% Obtain the fitted time from fit function
t50 = fzero(t50_time_evaluator, [0 1E6]);

end