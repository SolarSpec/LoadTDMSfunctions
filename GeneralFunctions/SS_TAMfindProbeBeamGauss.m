function [FitProbeBeam, ProbeBeamParameters, PinholePosition, FitDCmap, FitResidual, ResidualMap]...
    = SS_TAMfindProbeBeamGauss(InputDClevel,InputPinholeDiameter)
%SS_TAMfindProbeBeam Determine the probe beam parameters (assuming 2D Gaussian) that
% reproduce the DClevel map from a TAM acquisition

%% Define Globals to pass data that isn't fitted between functions
global PinholeDiameter ExpDClevel edgeLength Xcoord Ycoord

%% Function details

% Fit function FitCoefficients definitions:
% FitCoefficients(1): Amplitude of Gaussian
% FitCoefficients(2): Width of X direction Gaussian
% FitCoefficients(3): Width of Y direction Gaussian
% FitCoefficients(4): baseline offset
% FitCoefficients(5): X center position of pinhole
% FitCoefficients(6): Y center position of pinhole

% Calulate the global values for the other functions
ExpDClevel = InputDClevel;
edgeLength = size(ExpDClevel, 1);
PinholeDiameter = InputPinholeDiameter;

% Get X,Y coordinates of each points
[Xcoord,Ycoord] = meshgrid(1:edgeLength);

% Setup the fit
InitialGuess = [0.1, edgeLength/4, edgeLength/4, 0, edgeLength/2, edgeLength/2];
lb = [0, 0, 0, -max(max(ExpDClevel)), 0, 0];
ub = [inf, edgeLength, edgeLength, max(max(ExpDClevel)), 2*edgeLength, 2*edgeLength];
fitOptions = optimoptions("fmincon","PlotFcn","optimplotfval","Display","off",...
    "StepTolerance",1E-10,"TypicalX",[max(max(ExpDClevel)),edgeLength/2,edgeLength/2,1E-3,edgeLength/2,edgeLength/2]);

ResidualsFit = @(x)GenResiduals(x);

%Setup the problem
problem = createOptimProblem("fmincon",x0=InitialGuess,objective=ResidualsFit,lb=lb,ub=ub,options=fitOptions);

% Run the fit (local minimum)
% Old method without problem: [BestFitParam,FitResidual] = fmincon(ResidualsFit,InitialGuess,[],[],[],[],lb,ub,[],fitOptions);
% [BestFitParam,FitResidual] = fmincon(problem);

% Run the fit with Global search
gs = GlobalSearch('Display','iter');
[BestFitParam,FitResidual] = run(gs,problem);

% Calculate desired outputs
ProbeBeamParameters = [BestFitParam(1),edgeLength/2,BestFitParam(2), ...
    edgeLength/2, BestFitParam(3), BestFitParam(4)];

FitProbeBeam = Gauss2D(ProbeBeamParameters);

PinholePosition = [BestFitParam(5),BestFitParam(6)];

% Start iterating for each pixel position to generate the modeled TAM
% DClevel map
FitDCmap = zeros(edgeLength);

for posIndex = 1:1:numel(Xcoord)

    % Generate 2D Gaussian
    TempGauss = Gauss2D([BestFitParam(1),Xcoord(posIndex),BestFitParam(2),...
       Ycoord(posIndex), BestFitParam(3),BestFitParam(4)]);

    % Generate pinhole filter
    TempFilter = GenPinhole([BestFitParam(5), BestFitParam(6)]);

    % Determine light going through pinhole
    FitDCmap(Xcoord(posIndex),Ycoord(posIndex)) = sum(TempGauss .* TempFilter,'all');
end

ResidualMap = FitDCmap - InputDClevel;

end


function Residuals = GenResiduals(FitCoefficients)
% SS_GenResiduals calculates the sum of the root mean squared differences
% between the experimental DClevel TAM map and the modeled data.

% FitCoefficients definitions:
% FitCoefficients(1): Amplitude of Gaussian
% FitCoefficients(2): Width of X direction Gaussian
% FitCoefficients(3): Width of Y direction Gaussian
% FitCoefficients(4): baseline offset
% FitCoefficients(5): X center position of pinhole
% FitCoefficients(6): Y center position of pinhole

%% Define Globals to pass data that isn't fitted between functions
global ExpDClevel edgeLength Xcoord Ycoord

%% Function details

FitDClevel = zeros(edgeLength);

% Start iterating for each pixel position
for posIndex = 1:1:numel(Xcoord)

    % Generate 2D Gaussian
    TempGauss = Gauss2D([FitCoefficients(1),Xcoord(posIndex),FitCoefficients(2),...
        Ycoord(posIndex),FitCoefficients(3),FitCoefficients(4),FitCoefficients(5)]);

    % Generate pinhole filter
    TempFilter = GenPinhole([FitCoefficients(5), FitCoefficients(6)]);

    % Determine light going through pinhole
    FitDClevel(Xcoord(posIndex),Ycoord(posIndex)) = sum(TempGauss .* TempFilter,'all');
end

Residuals = sum(sqrt((FitDClevel - ExpDClevel).^2),'all');

end

function FilterOut = GenPinhole(PinholeCenter)
%GenPinhole Generate a 2D matrix of 0 and 1 based on a circular
%object at the center position.

%% Define Globals to pass data that isn't fitted between functions
global PinholeDiameter Xcoord Ycoord

%% Function details

centerX = PinholeCenter(1);
centerY = PinholeCenter(2);
PinholeRadius = PinholeDiameter/2;

% Determine which pixels are within the pinhole
FilterOut = (Xcoord - centerX).^2 + (Ycoord - centerY).^2 <= PinholeRadius.^2;

end


function GaussOut = Gauss2D(parameters)
%SS_2DGAUSS generates a 2D matrix formed from the vectorial multiplication
% of Gaussians in the X and Y directions. The edgeLength dictates the
% length of each side of the square.

%% Define Globals to pass data that isn't fitted between functions
global edgeLength

%% Function details
% Parameters definitions:
% parameters(1): Amplitude of Gaussian
% parameters(2): Center position of X direction Gaussian
% parameters(3): Width of X direction Gaussian
% parameters(4): Center position of Y direction Gaussian
% parameters(5): Width of Y direction Gaussian
% parameters(6): baseline offset

% Setup X and Y position vectors
Xvector = 1:1:edgeLength;
Yvector = 1:1:edgeLength;

% Generate X and Y direction vectors
XGauss = exp(-((Xvector-parameters(2))/parameters(3)).^2);
YGauss = transpose(exp(-((Yvector-parameters(4))/parameters(5)).^2));

% Generate 2D data from vectorial multiplication of X and Y gaussians, with
% an offset.
GaussOut = parameters(1) .* (YGauss * XGauss) + parameters(6);


end


