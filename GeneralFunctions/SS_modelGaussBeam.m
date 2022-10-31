function [Xfwhm,Yfwhm,avgFWHM,GaussBeamData] = SS_modelGaussBeam(parameters,GridSize,PixelSize)
%SS_MODELGAUSSBEAM Take the ProbeBeamParameters output of the SS_TAMfindProbeBeamGauss function.
% Calculates the FWHM (in length scale units) and generates a 3D dataset to plot the beam shape.

% Parameters definitions:
% parameters(1): Amplitude of Gaussian
% parameters(2): Center position of X direction Gaussian
% parameters(3): Width of X direction Gaussian, in pixels
% parameters(4): Center position of Y direction Gaussian
% parameters(5): Width of Y direction Gaussian, in pixels
% parameters(6): baseline offset

GaussBeamData = zeros(GridSize,GridSize,3);

% Setup X and Y position vectors
Xvector = linspace(0,2*parameters(2),GridSize);
Yvector = linspace(0,2*parameters(4),GridSize);

% Generate X and Y direction vectors
XGauss = exp(-((Xvector-parameters(2))/parameters(3)).^2);
YGauss = transpose(exp(-((Yvector-parameters(4))/parameters(5)).^2));

% Generate 2D data from vectorial multiplication of X and Y gaussians, with
% an offset.
GaussBeamData(:,:,3) = parameters(1) .* (YGauss * XGauss) + parameters(6);

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
