function [Fit2D, BestFitParameters, BestFitResiduals, BestFitDiffMap] = SS_2DFlatGaussFit(Input2Ddata)

% Setup X and Y position vectors
SizeInput = size(Input2Ddata);
Xvector = 1:1:SizeInput(1);
Yvector = 1:1:SizeInput(2);

% Find initial guesses
[MaxValue,MaxIndex] = max(Input2Ddata,[],'all');
[MaxCoordinatesY,MaxCoordinatesX] = ind2sub(size(Input2Ddata),MaxIndex);

InitialGuess = [20, 2, MaxCoordinatesX, ...
  1, 2, MaxCoordinatesY, 1, 0];

% Generate fit function so it can optimize the Gaussian parameters and use
% the input 2D data.
Fun4Fit = @(parameters)SS_2DFlatGauss(Input2Ddata,parameters);

% Tune fit options to help find the right fit.
FitLB = [0, 0, 0, 0, 0, 0, 0, 0];
FitUB = [inf, SizeInput(1), SizeInput(1), SizeInput(1), ...
    SizeInput(2), SizeInput(2), SizeInput(2), MaxValue];
FitOptions = optimset('PlotFcns',@optimplotfval);

% Perform the fit. Need to use fmincon for the Flat Top Gaussian to better
% sample the parameter space.
BestFitParameters = fmincon(Fun4Fit,InitialGuess,[], [], [], [], FitLB, FitUB, [], FitOptions);


% Generate X and Y direction vectors. Functional form taken from
% https://stats.stackexchange.com/questions/203629/is-there-a-plateau-shaped-distribution 
XFitFlatGauss = 1/(4*BestFitParameters(2))*(erf((Xvector-BestFitParameters(3)+BestFitParameters(2))/(BestFitParameters(4)*sqrt(2))) ...
    - erf((Xvector-BestFitParameters(3)-BestFitParameters(2))/(BestFitParameters(4)*sqrt(2))));
YFitFlatGauss = transpose(1/(4*BestFitParameters(5))*(erf((Yvector-BestFitParameters(6)+BestFitParameters(5))/(BestFitParameters(7)*sqrt(2))) ...
    - erf((Yvector-BestFitParameters(6)-BestFitParameters(5))/(BestFitParameters(7)*sqrt(2)))));

% Generate 2D data from vectorial multiplication of X and Y flat top gaussians, with
% an offset. Needs to be in the order of column vector * row vector to get
% appropriate outer matrix product.
Fit2D = BestFitParameters(1).*(YFitFlatGauss * XFitFlatGauss) + BestFitParameters(8);

% Calculate difference betewen modeled input and modeled data, then
% calculate the residuals.
BestFitDiffMap = Fit2D - Input2Ddata;
BestFitResiduals = sum(sum((BestFitDiffMap.^2)));


end

function [residuals] = SS_2DFlatGauss(Input2Ddata,parameters)
%SS_2DFlatGauss Calculate residuals from the square of the difference between a
% 2D matrix data input and the output from the vectorial multiplication of 
% a Gaussian in the X and Y directions.

% Parameters definitions:
% parameters(1): Amplitude of Flat Top Gaussian
% parameters(2): Width of the flat top in the X direction
% parameters(3): X direction center position
% parameters(4): Gaussian tail width of X direction Flat Top Gaussian 
% parameters(5): Width of the flat top in the Y direction 
% paremeters(6): Y direction center position
% paremeters(7): Gaussian tail width of Y direction Flat Top Gaussian
% parameters(8): baseline offset 

% Setup X and Y position vectors. Make bigger to allow for shift of flat
% top function which is centered to 0.
SizeInput = size(Input2Ddata);
Xvector = 1:1:SizeInput(1);
Yvector = 1:1:SizeInput(2);

% Generate X and Y direction vectors. Functional form taken from
% https://stats.stackexchange.com/questions/203629/is-there-a-plateau-shaped-distribution 
XFlatGauss = 1/(4*parameters(2))*(erf((Xvector-parameters(3)+parameters(2))/(parameters(4)*sqrt(2))) ...
    - erf((Xvector-parameters(3)-parameters(2))/(parameters(4)*sqrt(2))));
YFlatGauss = transpose(1/(4*parameters(5))*(erf((Yvector-parameters(6)+parameters(5))/(parameters(7)*sqrt(2))) ...
    - erf((Yvector-parameters(6)-parameters(5))/(parameters(7)*sqrt(2)))));

% Generate 2D data from vectorial multiplication of X and Y flat top gaussians, with
% an offset.
ModelFlatGauss = parameters(1).*(YFlatGauss * XFlatGauss) + parameters(8);

% Calculate difference betewen modeled input and modeled data, then
% calculate the residuals.
Diff = ModelFlatGauss - Input2Ddata;
residuals = sum(sum((Diff.^2)));

end

