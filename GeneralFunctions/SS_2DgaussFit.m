function [Fit2D, BestFitParameters, BestFitResiduals, BestFitDiffMap] = SS_2DgaussFit(Input2Ddata)

% Setup X and Y position vectors
SizeInput = size(Input2Ddata);
Xvector = 1:1:SizeInput(1);
Yvector = 1:1:SizeInput(2);

% Find initial guesses
[MaxValue,MaxIndex] = max(Input2Ddata,[],'all');
[MaxCoordinatesY,MaxCoordinatesX] = ind2sub(size(Input2Ddata),MaxIndex);

InitialGuess = [sqrt(MaxValue), MaxCoordinatesX, 1, sqrt(MaxValue), ...
    MaxCoordinatesY, 1, 0];

% Generate fit function so it can optimize the Gaussian parameters and use
% the input 2D data.
Fun4Fit = @(parameters)SS_2Dgauss(Input2Ddata,parameters);

% Perform the fit. 
BestFitParameters = fminsearch(Fun4Fit,InitialGuess);

% Generate X and Y direction vectors
XFitGauss = BestFitParameters(1)*exp(-((Xvector-BestFitParameters(2))/BestFitParameters(3)).^2);
YFitGauss = transpose(BestFitParameters(4)*exp(-((Yvector-BestFitParameters(5))/BestFitParameters(6)).^2));

% Generate 2D data from vectorial multiplication of X and Y gaussians, with
% an offset. Needs to be in the order of column vector * row vector to get
% appropriate outer matrix product.
Fit2D = YFitGauss * XFitGauss + BestFitParameters(7);

% Calculate difference betewen modeled input and modeled data, then
% calculate the residuals.
BestFitDiffMap = Fit2D - Input2Ddata;
BestFitResiduals = sum(sum((BestFitDiffMap.^2)));


end

function [residuals] = SS_2Dgauss(Input2Ddata,parameters)
%SS_2DGAUSS Calculate residuals from the square of the difference between a
% 2D matrix data input and the output from the vectorial multiplication of 
% a Gaussian in the X and Y directions.

% Parameters definitions:
% parameters(1): Amplitude of X direction Gaussian
% parameters(2): Center position of X direction Gaussian
% parameters(3): Width of X direction Gaussian
% parameters(4): Amplitude of Y direction Gaussian
% parameters(5): Center position of Y direction Gaussian
% parameters(6): Width of Y direction Gaussian
% parameters(7): baseline offset 

% Setup X and Y position vectors
SizeInput = size(Input2Ddata);
Xvector = 1:1:SizeInput(1);
Yvector = 1:1:SizeInput(2);

% Generate X and Y direction vectors
XGauss = parameters(1)*exp(-((Xvector-parameters(2))/parameters(3)).^2);
YGauss = transpose(parameters(4)*exp(-((Yvector-parameters(5))/parameters(6)).^2));

% Generate 2D data from vectorial multiplication of X and Y gaussians, with
% an offset.
ModelGauss = YGauss * XGauss + parameters(7);

% Calculate difference betewen modeled input and modeled data, then
% calculate the residuals.
Diff = ModelGauss - Input2Ddata;
residuals = sum(sum((Diff.^2)));

end

