%% Input description
% The workspace should contain the following:
% A row vector of Data, which is your TRPL decay data
% A row vector of Time, which serves as the X axis info. The units are assumed to be in nanoseconds
% A row vector of IRF, which is the TRPL decay of the IRF.

% The X axis for Data and IRF need to be the same. All vectors need to be of the same length. Only fit one trace at a time.

%% Parameters that the user might want to change

% Toggle figure export on or off
FigExport = false;

% Define fit window
WindowLB = 48;
WindowUB = 105;

% Define area to get baseline level guess
WindowBaseline = [900:1100];

%% Choose the fit model
FitList = {'1 power law', '1 exponential, 1 power law','2 exponentials, 1 power law',...
    '1 exponential, 1 second order','2 exponentials, 1 second order'};
[ListChoice,tf] = listdlg('ListString',FitList,'SelectionMode','single',...
    'ListSize',[200,75],'Name','Choose fitting model');

if tf==0
    error("Didn't choose a fit model. Exiting...")
end

switch ListChoice
    case 1
        FitModel = '1power';
    case 2
        FitModel = '1exp1power';
    case 3
        FitModel = '2exp1power';
    case 4
        FitModel = '1exp1second';
    case 5
        FitModel = '2exp1second';
end

disp(['Fit model chosen: ', FitModel])

%% Start the parallel pool
gcp;


%% Get rid of the trailling 0 points
while Data(end) == 0
    Data(end) = [];
    IRF(end) = [];
    Time(end) = [];
end

%% Prepare for fit data at the different IRF shifts

switch FitModel
    case '1power'       % (1 exp, 1 power law parameters)
        numParam = 5;
    case '1exp1power'   % (1 exp, 1 power law parameters)
        numParam = 7;
    case '2exp1power'   % (2 exp, 1 power law parameters)
        numParam = 8;
    case '1exp1second'  % (1 exp, 1 second order components)
        numParam = 6;
    case '2exp1second'  % (2 exp, 1 second order components)
        numParam = 8;
end

%% Restrict to fit window

TimeZero = Time(max(Data)==Data);

WindowData = Data(Time>WindowLB & Time<WindowUB);
WindowIRF = IRF(Time>WindowLB & Time<WindowUB);
WindowTime = Time(Time>WindowLB & Time<WindowUB);

%% Setup initial parameters

switch FitModel
    case '1power'       % (1 power law parameters)

        %param(1);      %Y offset
        %param(2);      %Amplitude of the Power law component
        %param(3);      %Power law component onset parameter
        %param(4);      %Power law exponent parameter (alpha)
        %param(5);      %IRF time shift

        InitialGuess = [mean(Data(WindowBaseline)); max(Data)/2; 0.1; 0.5; 0];
        LowerBound = [0; 0; 0; 0; -3];
        UpperBound = [mean(Data(WindowBaseline))*5; inf; inf; 5; 3];

        TRPLfitfun = @(param)Fit_IRF_1power(param,WindowTime,WindowData,WindowIRF,TimeZero);

    case '1exp1power'       % (1 exp, 1 power law parameters)

        %param(1);      %Y offset
        %param(2);      %Amplitude of the exponential component
        %param(3);      %Exponential component lifetime
        %param(4);      %Amplitude of the Power law component
        %param(5);      %Power law component onset parameter
        %param(6);      %Power law exponent parameter (alpha)
                %param(7);      %IRF time shift

        InitialGuess = [mean(Data(WindowBaseline)); max(Data)/2; 2 ; max(Data)/2; 0.765; 1.75; 0];
        LowerBound = [0; 0; 0; 0; 0.765; 1.75; -3];
        UpperBound = [mean(Data(WindowBaseline))*5; inf; 2*max(Time) ; inf; 0.765; 1.75; 3];

        TRPLfitfun = @(param)Fit_IRF_1exp1power(param,WindowTime,WindowData,WindowIRF,TimeZero);

    case '2exp1power'    % (2 exp, 1 power law parameters)

        %param(1);      %Y offset
        %param(2);      %Amplitude of the first exponential component
        %param(3);      %First exponential component lifetime
        %param(4);      %Amplitude of the second exponential component
        %param(5);      %Second exponential component lifetime
        %param(6);      %Amplitude of the Power law component
        %param(7);      %Power law component onset parameter
        %param(8);      %Power law exponent parameter (alpha)
                %param(9);      %IRF time shift

        InitialGuess = [mean(Data(WindowBaseline)); max(Data)/2; 0.5 ; max(Data)/2; 2; max(Data); 0.88; 0.5; 0];
        LowerBound = [0; 0; 0; 0; 0; 0; 0.86; 1.74; -3];
        UpperBound = [mean(Data(WindowBaseline))*5; inf; 0.01; inf; 2*max(Time); inf; 0.89; 1.76; 3];

        TRPLfitfun = @(param)Fit_IRF_2exp1power(param,WindowTime,WindowData,WindowIRF,TimeZero);

    case '1exp1second'  % (1 exp, 1 second order components)
        %param(1);         %Y offset
        %param(2);         %Amplitude of the first exponential component
        %param(3);         %First exponential component lifetime
        %param(4);         %Amplitude of the second order component
        %param(5);         %Second order kinetic parameter
                %param(6);      %IRF time shift

        InitialGuess = [mean(Data(WindowBaseline)); max(Data)/2; 2 ; max(Data); 0.1; 0];
        LowerBound = [0; 0; 0.2; 0; 0.001; -3];
        UpperBound = [mean(Data(WindowBaseline))*5; inf; 2*max(Time); inf; inf; 3];

        TRPLfitfun = @(param)Fit_IRF_1exp1second(param,WindowTime,WindowData,WindowIRF,TimeZero);

    case '2exp1second'    % (2 exp, 1 power law parameters)

        %param(1);      %Y offset
        %param(2);      %Amplitude of the first exponential component
        %param(3);      %First exponential component lifetime
        %param(4);      %Amplitude of the second exponential component
        %param(5);      %Second exponential component lifetime
        %param(6);      %Amplitude of the second order component
        %param(7);      %Second order kinetic parameter
                %param(8);      %IRF time shift

        InitialGuess = [mean(Data(WindowBaseline)); max(Data)/2; 0.4 ; max(Data)/2; 0.8; max(Data); 0.01; 0];
        LowerBound = [0; 0; 0.2; 0; 0.2; 0; 0.001; -3];
        UpperBound = [mean(Data(WindowBaseline))*5; inf; 8; inf; 8; inf; 0.03; 3];

        TRPLfitfun = @(param)Fit_IRF_2exp1second(param,WindowTime,WindowData,WindowIRF,TimeZero);
end

%% Global minimization

options = optimoptions(@fmincon);
problem = createOptimProblem('fmincon','objective',TRPLfitfun,'x0',InitialGuess,...
    'lb',LowerBound,'ub',UpperBound);
problem.options.Display = 'none';
% %problem.options.PlotFcns = @optimplotfval;
% %problem.options.PlotFcns = @optimplotx;

%[x,fval,eflag,output] = fmincon(problem); %Local minimization for testing problem structure.

ms = MultiStart;
ms.UseParallel = true;
ms.StartPointsToRun = 'all';
ms.XTolerance = 1E-12;
ms.FunctionTolerance = 1E-12;
numMultiStartRuns = 50;

%gs = GlobalSearch(ms);
%rng(14,'twister')   % for reproducibility
tic                 % for timing
%[paramgs,fvalgs] = run(gs,problem);
[paramgs,fvalgs] = run(ms,problem,numMultiStartRuns);
toc


%% Generate best fit data

switch FitModel
    case '1power'
        [FinalFit, FitYdataPower] = Decay_IRF_1power(paramgs,Time,IRF,TimeZero);
    case '1exp1power'
        [FinalFit, FitYdataExp, FitYdataPower] = Decay_IRF_1exp1power(paramgs,Time,IRF,TimeZero);
    case '2exp1power'
        [FinalFit, FitYdataExp, FitYdataExp2, FitYdataPower] = Decay_IRF_2exp1power(paramgs,Time,IRF,TimeZero);
    case '1exp1second'
        [FinalFit, FitYdataExp, FitYdataSecond] = Decay_IRF_1exp1second(paramgs,Time,IRF,TimeZero);
    case '2exp1second'
        [FinalFit, FitYdataExp, FitYdataExp2, FitYdataSecond] = Decay_IRF_2exp1second(paramgs,Time,IRF,TimeZero);
end

CleanWindowData = WindowData;
CleanWindowData(WindowData==0) = 1;         % Remove 0's as can't divide by 0 to calculate the reduced Chi-squared.

CleanData = Data;
CleanData(Data==0) = 1;                     % Remove 0's as can't divide by 0 to calculate the reduced Chi-squared.

FitDelta = CleanData - FinalFit;
FitRes = FitDelta./sqrt(CleanData);
ChiSq = sum(FitRes(Time>WindowLB & Time<WindowUB).^2./CleanWindowData);
RedChiSq = ChiSq/(length(CleanWindowData)-length(InitialGuess)-1);       % Length Initial Guess gives you the number of fitting parameters.


%% Plot ouputs
subplot(10,1,1:6)

semilogy(Time,Data)
hold on
semilogy(Time,FinalFit,'LineWidth',0.75)
yline(paramgs(1),'--k')
ylabel('Counts')
fontsize(gca,12,"points")
switch FitModel
    case '1power'
        plot(Time,FitYdataPower)
        legend({'Data','Total Fit','Baseline','Power law comp.'},...
            'Box','off','FontSize',8,'NumColumns',2)

    case '1exp1power'
        plot(Time,FitYdataExp)
        plot(Time,FitYdataPower)
        legend({'Data','Total Fit','Baseline','Exponential comp.',...
            'Power law comp.'},'Box','off','FontSize',8,'NumColumns',2)

    case '2exp1power'
        plot(Time,FitYdataExp)
        plot(Time,FitYdataExp2)
        plot(Time,FitYdataPower)
        legend({'Data','Total Fit','Baseline','Exponential comp. 1',...
            'Exponential comp. 2','Power law comp.'},'Box','off','FontSize',8,'NumColumns',2)

    case '1exp1second'
        plot(Time,FitYdataExp)
        plot(Time,FitYdataSecond)
        legend({'Data','Total Fit','Baseline','Exponential comp.',...
            'Second ord. comp.'},'Box','off','FontSize',8,'NumColumns',2)

    case '2exp1second'
        plot(Time,FitYdataExp)
        semilogy(Time,FitYdataExp2)
        semilogy(Time,FitYdataSecond)
        legend({'Data','Total Fit','Baseline','Exponential comp. 1',...
            'Exponential comp. 2','Second ord. comp.'},'Box','off','FontSize',8,'NumColumns',2)

end


hold off
subplot(10,1,8:10)
plot(Time,FitRes)
yline(0)
xlabel('Time (ns)')
ylabel('Residuals')
fontsize(gca,12,"points")

if FigExport == false
    title(['IRF shift of ' num2str(paramgs(end),3), '. Chi-squared = ' num2str(ChiSq,3), '. Reduced Chi-squared = ' num2str(RedChiSq,3)])
else
    [Figfile,Figpath] = uiputfile('.emf','Export figure to...');
    exportgraphics(gcf,fullfile(Figpath,Figfile))
end

%% Output additional parameters from the fit

% Find t50% from fit
SearchTime = Time(Time>TimeZero & Time<WindowUB);
SearchData = FinalFit(Time>TimeZero & Time<WindowUB);

HalfIntTimes = SearchTime(min(abs(SearchData-max(SearchData)/2))==abs(SearchData-max(SearchData)/2));
RoughHalfIntTime = mean(HalfIntTimes);          %In case there are multiple values

DataSpacing = Time(end)-Time(end-1);
FineTime = RoughHalfIntTime-2*DataSpacing:DataSpacing/50:RoughHalfIntTime+2*DataSpacing;
FineData = interp1(Time,Data,FineTime);
FinalHalfIntTime = FineTime(min(abs(FineData-max(SearchData)/2))==abs(FineData-max(SearchData)/2));

t50 = FinalHalfIntTime(1) - TimeZero;              %TimeZero is defined as the time where the Ydata is max

% Find relative emission from the different components

switch FitModel
    case '1power'

        RelEmPower = 1;         % Trivial since there is only 1 component.

    case '1exp1power'

        ExpCompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataExp(Time>WindowLB & Time<WindowUB));
        PowerCompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataPower(Time>WindowLB & Time<WindowUB));

        RelEmExp = ExpCompIntegral/(ExpCompIntegral+PowerCompIntegral);
        RelEmPower = PowerCompIntegral/(ExpCompIntegral+PowerCompIntegral);

    case '2exp1power'

        Exp1CompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataExp(Time>WindowLB & Time<WindowUB));
        Exp2CompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataExp2(Time>WindowLB & Time<WindowUB));
        PowerCompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataPower(Time>WindowLB & Time<WindowUB));

        RelEmExp1 = Exp1CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);
        RelEmExp2 = Exp2CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);
        RelEmPower = PowerCompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);

    case '1exp1second'

        ExpCompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataExp(Time>WindowLB & Time<WindowUB));
        SecondCompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataSecond(Time>WindowLB & Time<WindowUB));

        RelEmExp = ExpCompIntegral/(ExpCompIntegral+SecondCompIntegral);
        RelEmSecond = SecondCompIntegral/(ExpCompIntegral+SecondCompIntegral);

    case '2exp1second'

        Exp1CompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataExp(Time>WindowLB & Time<WindowUB));
        Exp2CompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataExp2(Time>WindowLB & Time<WindowUB));
        SecondCompIntegral = trapz(Time(Time>WindowLB & Time<WindowUB),FitYdataSecond(Time>WindowLB & Time<WindowUB));

        RelEmExp1 = Exp1CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);
        RelEmExp2 = Exp2CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);
        RelEmSecond = SecondCompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);

end


%% Setup fit function (1 power law)

function RedChiSq = Fit_IRF_1power(param,Time,Data,IRF,TimeZero)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
PowerAmp = param(2);        %Amplitude of the Power law component
PowerOnset = param(3);      %Power law component onset parameter
PowerAlpha = param(4);      %Power law exponent parameter (alpha)
IRFshift = param(5);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");


% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
YdataPL = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

GuessKin = Y0 + YdataPL;

nonZeroInd = find(Data);
Delta = (Data(nonZeroInd) - GuessKin(nonZeroInd));  %Can't include where Data = 0 because of the division
ChiSq = sum((Delta.^2./Data(nonZeroInd)));          %Works much better than least squares minimization Res = sum(sum( (Data - GuessKin).^2));
RedChiSq = ChiSq/(length(Data)-length(param)-1);    %Compare reduced ChiSq to better compare different datasets and fit windows

end

%% Setup fit function (1 exp, 1 power law)

function RedChiSq = Fit_IRF_1exp1power(param,Time,Data,IRF,TimeZero)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
PowerAmp = param(4);        %Amplitude of the Power law component
PowerOnset = param(5);      %Power law component onset parameter
PowerAlpha = param(6);      %Power law exponent parameter (alpha)
IRFshift = param(7);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");


% Generate Data for exponential component
TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);


% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
YdataPL = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

GuessKin = Y0 + YdataExp + YdataPL;

nonZeroInd = find(Data);
Delta = (Data(nonZeroInd) - GuessKin(nonZeroInd));  %Can't include where Data = 0 because of the division
ChiSq = sum((Delta.^2./Data(nonZeroInd)));          %Works much better than least squares minimization Res = sum(sum( (Data - GuessKin).^2));
RedChiSq = ChiSq/(length(Data)-length(param)-1);    %Compare reduced ChiSq to better compare different datasets and fit windows

end

%% Setup fit function (2 exp, 1 power law)

function RedChiSq = Fit_IRF_2exp1power(param,Time,Data,IRF,TimeZero)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
PowerAmp = param(6);        %Amplitude of the Power law component
PowerOnset = param(7);      %Power law component onset parameter
PowerAlpha = param(8);      %Power law exponent parameter (alpha)
IRFshift = param(9);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
YdataExp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
YdataPL = ConvYPL(1,NnegX+1:length(t)+NnegX);

% Calculate output

GuessKin = Y0 + YdataExp + YdataExp2 + YdataPL;

nonZeroInd = find(Data);
Delta = (Data(nonZeroInd) - GuessKin(nonZeroInd));  %Can't include where Data = 0 because of the division
ChiSq = sum((Delta.^2./Data(nonZeroInd)));          %Works much better than least squares minimization Res = sum(sum( (Data - GuessKin).^2));
RedChiSq = ChiSq/(length(Data)-length(param)-1);    %Compare reduced ChiSq to better compare different datasets and fit windows

end

%% Setup fit function (1 exp, 1 second order)

function RedChiSq = Fit_IRF_1exp1second(param,Time,Data,IRF,TimeZero)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
SecondAmp = param(4);       %Amplitude of the second order component
SecondKin = param(5);       %Second order kinetic parameter
IRFshift = param(6);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for Second Order Component
TempYSecond = SecondAmp./((1+(t-TimeZero)/SecondKin));
TempYSecond(t<TimeZero) = 0;
ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
YdataSecond = ConvYSecond(1,NnegX+1:length(t)+NnegX);

% Calculate output

GuessKin = Y0 + YdataExp + YdataSecond;

nonZeroInd = find(Data);
Delta = (Data(nonZeroInd) - GuessKin(nonZeroInd));  %Can't include where Data = 0 because of the division
ChiSq = sum((Delta.^2./Data(nonZeroInd)));          %Works much better than least squares minimization Res = sum(sum( (Data - GuessKin).^2));
RedChiSq = ChiSq/(length(Data)-length(param)-1);    %Compare reduced ChiSq to better compare different datasets and fit windows
end

%% Setup fit function (2 exp, 1 second order)

function RedChiSq = Fit_IRF_2exp1second(param,Time,Data,IRF,TimeZero)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
SecondAmp = param(6);       %Amplitude of the second order component
SecondKin = param(7);       %Second order kinetic parameter
IRFshift = param(8);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
YdataExp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
YdataExp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Generate Data for Second Order Component
TempYSecond = SecondAmp./((1+(t-TimeZero)/SecondKin));
TempYSecond(t<TimeZero) = 0;
ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
YdataSecond = ConvYSecond(1,NnegX+1:length(t)+NnegX);


% Calculate output

GuessKin = Y0 + YdataExp + YdataExp2 + YdataSecond;

nonZeroInd = find(Data);
Delta = (Data(nonZeroInd) - GuessKin(nonZeroInd));  %Can't include where Data = 0 because of the division
ChiSq = sum((Delta.^2./Data(nonZeroInd)));          %Works much better than least squares minimization Res = sum(sum( (Data - GuessKin).^2));
RedChiSq = ChiSq/(length(Data)-length(param)-1);    %Compare reduced ChiSq to better compare different datasets and fit windows

end

%% Setup decay function for getting the convolved fitted data (1 power)

function [TotalOutput,PowerComp] = Decay_IRF_1power(param,Time,IRF,TimeZero)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
PowerAmp = param(2);        %Amplitude of the Power law component
PowerOnset = param(3);      %Power law component onset parameter
PowerAlpha = param(4);      %Power law exponent parameter (alpha)
IRFshift = param(5);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");


% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
PowerComp = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + PowerComp;

end

%% Setup decay function for getting the convolved fitted data (1 exp, 1 power)

function [TotalOutput,ExpComp,PowerComp] = Decay_IRF_1exp1power(param,Time,IRF,TimeZero)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
PowerAmp = param(4);        %Amplitude of the Power law component
PowerOnset = param(5);      %Power law component onset parameter
PowerAlpha = param(6);      %Power law exponent parameter (alpha)
IRFshift = param(7);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for exponential component
TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);


% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
PowerComp = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + ExpComp + PowerComp;

end

%% Setup decay function for getting the convolved fitted data (2 exp, 1 power)

function [TotalOutput,ExpComp,ExpComp2,PowerComp] = Decay_IRF_2exp1power(param,Time,IRF,TimeZero)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
PowerAmp = param(6);        %Amplitude of the Power law component
PowerOnset = param(7);      %Power law component onset parameter
PowerAlpha = param(8);      %Power law exponent parameter (alpha)
IRFshift = param(9);        %Time shift of IRF

%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
ExpComp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Generate Data for Power Law Component
TempYPL = PowerAmp./(((1+(t-TimeZero)/PowerOnset)).^PowerAlpha);
TempYPL(t<TimeZero) = 0;
ConvYPL = conv(TempYPL,ShiftedIRF/sum(ShiftedIRF));
PowerComp = ConvYPL(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + ExpComp + ExpComp2 + PowerComp;

end

%% Setup decay function for getting the convolved fitted data (1 exp, 1 second)

function [TotalOutput,ExpComp,SecondComp] = Decay_IRF_1exp1second(param,Time,IRF,TimeZero)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
SecondAmp = param(4);       %Amplitude of the second ordercomponent
SecondKin = param(5);       %Second order component kinetic parameter
IRFshift = param(6);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');

%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for exponential component
TempYExp = ExpAmp.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);


% Generate Data for second order component
TempYSecond = SecondAmp./(1+(t-TimeZero)/SecondKin);
TempYSecond(t<TimeZero) = 0;
ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
SecondComp = ConvYSecond(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + ExpComp + SecondComp;

end

%% Setup decay function for getting the convolved fitted data (2 exp, 1 second)

function [TotalOutput,ExpComp,ExpComp2,SecondComp] = Decay_IRF_2exp1second(param,Time,IRF,TimeZero)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
SecondAmp = param(6);       %Amplitude of the second ordercomponent
SecondKin = param(7);       %Second order component kinetic parameter
IRFshift = param(8);        %Time shift of IRF


%Find the index of the first component of where the 'non-zero' time data is
NnegX = find(t<TimeZero,1,'last');


%Get the IRF values at the shifted time. Use "pchip" to also extrapolate
%when outside the range
ShiftedIRF = interp1(t,IRF,t+IRFshift,"pchip");

% Generate Data for the first exponential component
TempYExp = ExpAmp1.*exp(-(t-TimeZero)./ExpTau1);
TempYExp(t<TimeZero) = 0;
ConvYExp = conv(TempYExp,ShiftedIRF/sum(ShiftedIRF));
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
ExpComp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Generate Data for second order component
TempYSecond = SecondAmp./(1+(t-TimeZero)/SecondKin);
TempYSecond(t<TimeZero) = 0;
ConvYSecond = conv(TempYSecond,ShiftedIRF/sum(ShiftedIRF));
SecondComp = ConvYSecond(1,NnegX+1:length(t)+NnegX);


% Calculate output

TotalOutput = Y0 + ExpComp + ExpComp2 + SecondComp;

end