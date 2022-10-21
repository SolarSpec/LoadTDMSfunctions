function [FitProbeBeam, ProbeBeamParameters, PinholePosition, FitDCmap, FitResidual, ResidualMap] ...
    = SS_TAMfindProbeBeamFlat(InputDClevel,InputPinholeDiameter)
%SS_TAMfindProbeBeam Determine the probe beam parameters (assuming 2D Gaussian) that
% reproduce the DClevel map from a TAM acquisition

%% Define Globals to pass data that isn't fitted between functions
global PinholeDiameter ExpDClevel edgeLength Xcoord Ycoord

%% Function details

% Fit function FitCoefficients definitions:
% FitCoefficients(1): Amplitude of Flat Top Gaussian
% FitCoefficients(2): Width of the flat top in the X direction
% FitCoefficients(3): Gaussian tail width of X direction Flat Top Gaussian 
% FitCoefficients(4): Width of the flat top in the Y direction 
% FitCoefficients(5): Gaussian tail width of Y direction Flat Top Gaussian
% FitCoefficients(6): baseline offset 
% FitCoefficients(7): X center position of pinhole
% FitCoefficients(8): Y center position of pinhole

% Calulate the global values for the other functions
ExpDClevel = InputDClevel;
edgeLength = size(ExpDClevel, 1);
PinholeDiameter = InputPinholeDiameter;

% Get X,Y coordinates of each points
[Xcoord,Ycoord] = meshgrid(1:edgeLength);

% Setup the fit
InitialGuess = [0.1, edgeLength/4, edgeLength/4, edgeLength/4, edgeLength/4,...
    0, edgeLength/2,edgeLength/2];
lb = [0, 0, 0, 0, 0, -max(max(ExpDClevel)), 0, 0];
ub = [inf, edgeLength, edgeLength, edgeLength, edgeLength, max(max(ExpDClevel)), 2*edgeLength, 2*edgeLength];
fitOptions = optimoptions("fmincon","PlotFcn","optimplotfval","Display","off",...
    "StepTolerance",1E-10,"TypicalX",...
    [max(max(ExpDClevel)),edgeLength/2,edgeLength/2,edgeLength/2,edgeLength/2,1E-3,edgeLength/2,edgeLength/2]);

ResidualsFit = @(x)GenResiduals(x);

%Setup the problem
problem = createOptimProblem("fmincon",x0=InitialGuess,objective=ResidualsFit,lb=lb,ub=ub,options=fitOptions);

% Run the fit (local minimum)
% Old method without problem: BestFitParam = fmincon(ResidualsFit,InitialGuess,[],[],[],[],lb,ub,[],fitOptions);
% [BestFitParam,~] = fmincon(problem);

% Run the fit with Global search
gs = GlobalSearch('Display','iter');
[BestFitParam,FitResidual] = run(gs,problem);

% Calculate desired outputs
ProbeBeamParameters = [BestFitParam(1),BestFitParam(2),edgeLength/2, ...
    BestFitParam(3), BestFitParam(4),edgeLength/2, BestFitParam(5), BestFitParam(6)];

FitProbeBeam = Flat2D(ProbeBeamParameters);

PinholePosition = [BestFitParam(7),BestFitParam(8)];

% Start iterating for each pixel position to generate the modeled TAM
% DClevel map
FitDCmap = zeros(edgeLength);

for posIndex = 1:1:numel(Xcoord)

    % Generate 2D Gaussian
    TempGauss = Flat2D([BestFitParam(1),BestFitParam(2),Xcoord(posIndex),...
               BestFitParam(3), BestFitParam(4),Ycoord(posIndex), BestFitParam(5),BestFitParam(6)]);

    % Generate pinhole filter
    TempFilter = GenPinhole([BestFitParam(7), BestFitParam(8)]);

    % Determine light going through pinhole
    FitDCmap(Xcoord(posIndex),Ycoord(posIndex)) = sum(TempGauss .* TempFilter,'all');
end

ResidualMap = FitDCmap - InputDClevel;

end


function Residuals = GenResiduals(FitCoefficients)
% SS_GenResiduals calculates the sum of the root mean squared differences
% between the experimental DClevel TAM map and the modeled data.

% FitCoefficients(1): Amplitude of Flat Top Gaussian
% FitCoefficients(2): Width of the flat top in the X direction
% FitCoefficients(3): Gaussian tail width of X direction Flat Top Gaussian 
% FitCoefficients(4): Width of the flat top in the Y direction 
% FitCoefficients(5): Gaussian tail width of Y direction Flat Top Gaussian
% FitCoefficients(6): baseline offset 
% FitCoefficients(7): X center position of pinhole
% FitCoefficients(8): Y center position of pinhole

%% Define Globals to pass data that isn't fitted between functions
global ExpDClevel edgeLength Xcoord Ycoord

%% Function details

FitDClevel = zeros(edgeLength);

% Start iterating for each pixel position
for posIndex = 1:1:numel(Xcoord)

    % Generate 2D Flat top Gaussian
    TempGauss = Flat2D([FitCoefficients(1),FitCoefficients(2), Xcoord(posIndex),...
        FitCoefficients(3),FitCoefficients(4),Ycoord(posIndex),...
        FitCoefficients(5),FitCoefficients(6)]);

    % Generate pinhole filter
    TempFilter = GenPinhole([FitCoefficients(7), FitCoefficients(8)]);

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


function FlatOut = Flat2D(parameters)
%SS_2DGAUSS generates a 2D matrix formed from the vectorial multiplication
% of Gaussians in the X and Y directions. The edgeLength dictates the
% length of each side of the square.

%% Define Globals to pass data that isn't fitted between functions
global edgeLength

%% Function details
% Parameters definitions:
% parameters(1): Amplitude of Flat Top Gaussian
% parameters(2): Width of the flat top in the X direction
% parameters(3): X direction center position
% parameters(4): Gaussian tail width of X direction Flat Top Gaussian 
% parameters(5): Width of the flat top in the Y direction 
% paremeters(6): Y direction center position
% paremeters(7): Gaussian tail width of Y direction Flat Top Gaussian
% parameters(8): baseline offset 

% Setup X and Y position vectors
Xvector = 1:1:edgeLength;
Yvector = 1:1:edgeLength;

% Generate X and Y direction vectors. Functional form taken from
% https://stats.stackexchange.com/questions/203629/is-there-a-plateau-shaped-distribution 
XFlatGauss = 1/(4*parameters(2))*(erf((Xvector-parameters(3)+parameters(2))/(parameters(4)*sqrt(2))) ...
    - erf((Xvector-parameters(3)-parameters(2))/(parameters(4)*sqrt(2))));
YFlatGauss = transpose(1/(4*parameters(5))*(erf((Yvector-parameters(6)+parameters(5))/(parameters(7)*sqrt(2))) ...
    - erf((Yvector-parameters(6)-parameters(5))/(parameters(7)*sqrt(2)))));

% Generate 2D data from vectorial multiplication of X and Y flat top gaussians, with
% an offset.
FlatOut = parameters(1).*(YFlatGauss * XFlatGauss) + parameters(8);


end


