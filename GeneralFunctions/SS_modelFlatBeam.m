function [Xfwhm,Yfwhm,avgFWHM,GaussBeamData] = SS_modelFlatBeam(parameters,GridSize,PixelSize)
%SS_MODELFlatBEAM Take the ProbeBeamParameters output of the SS_TAMfindProbeBeamFlats function.
% Calculates the FWHM (in length scale units) and generates a 3D dataset to plot the beam shape.

% Parameters definitions:
% parameters(1): Amplitude of Flat Top Gaussian
% parameters(2): Width of the flat top in the X direction, in pixels
% parameters(3): X direction center position
% parameters(4): Gaussian tail width of X direction Flat Top Gaussian, in pixels
% parameters(5): Width of the flat top in the Y direction, in pixels 
% paremeters(6): Y direction center position
% paremeters(7): Gaussian tail width of Y direction Flat Top Gaussian, in pixels
% parameters(8): baseline offset 

GaussBeamData = zeros(GridSize,GridSize,3);

% Setup X and Y position vectors
Xvector = linspace(0,2*parameters(3),GridSize);
Yvector = linspace(0,2*parameters(6),GridSize);

% Generate X and Y direction vectors. Functional form taken from
% https://stats.stackexchange.com/questions/203629/is-there-a-plateau-shaped-distribution 
XFlatGauss = 1/(4*parameters(2))*(erf((Xvector-parameters(3)+parameters(2))/(parameters(4)*sqrt(2))) ...
    - erf((Xvector-parameters(3)-parameters(2))/(parameters(4)*sqrt(2))));
YFlatGauss = transpose(1/(4*parameters(5))*(erf((Yvector-parameters(6)+parameters(5))/(parameters(7)*sqrt(2))) ...
    - erf((Yvector-parameters(6)-parameters(5))/(parameters(7)*sqrt(2)))));

% Generate 2D data from vectorial multiplication of X and Y flat top gaussians, with
% an offset.
GaussBeamData(:,:,3) = parameters(1).*(YFlatGauss * XFlatGauss) + parameters(8);

% Generate mesh grid and convert to length units
[MeshX,MeshY] = meshgrid(Xvector,Yvector);

GaussBeamData(:,:,1) = MeshX*PixelSize;
GaussBeamData(:,:,2) = MeshY*PixelSize;

% Find the max to determine the FWHM from
BeamAmp = GaussBeamData(:,:,3);
[~,MaxIndex] = max(BeamAmp(:));
[Xindex,Yindex] = ind2sub(size(GaussBeamData(:,:,3)),MaxIndex);

% Calculate the FWHM in the X and Y directions and take the average.
[~,~,Xfwhm,~] = findpeaks(GaussBeamData(Xindex,:,3),GaussBeamData(Xindex,:,1));
[~,~,Yfwhm,~] = findpeaks(GaussBeamData(:,Yindex,3),GaussBeamData(:,Yindex,2));

avgFWHM = (Xfwhm+Yfwhm)/2;

% Plot output
surf(GaussBeamData(:,:,1),GaussBeamData(:,:,2),GaussBeamData(:,:,3));
shading interp
end
