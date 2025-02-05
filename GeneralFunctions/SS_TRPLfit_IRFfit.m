%% Input description
% The workspace should contain the following:
% A row vector of RawData, which is your TRPL decay data
% A row vector of RawTime, which serves as the X axis info. The units are assumed to be in nanoseconds
% A row vector of RawIRF, which is the TRPL decay of the IRF.

% The X axis for Data and IRF need to be the same. All vectors need to be of the same length. Only fit one trace at a time.

Time = RawTime; % Define from raw to avoid issues of different vector lengths when just changing one of them and there are trailling zeros.
IRF = RawIRF;
Data = RawData;

%% Parameters that the user might want to change

% Toggle figure export on or off
FigExport = false;

% Define fit window. %Editable.
WindowLB = 48;
WindowUB = 80;

% Define area to get baseline level guess. %Editable.
WindowBaseline = 900:1100;

%% Choose the fit model
FitList = {'1 exponential', '1 power law', '1 exponential, 1 power law',...
    '2 exponentials','2 exponentials, 1 power law',...
    '1 exponential, 1 second order','2 exponentials, 1 second order'};
[ListChoice,tf] = listdlg('ListString',FitList,'SelectionMode','single',...
    'ListSize',[200,100],'Name','Choose fitting model');

if tf==0
    error("Didn't choose a fit model. Exiting...")
end

switch ListChoice
    case 1
        FitModel = '1exp';
    case 2
        FitModel = '1power';
    case 3
        FitModel = '1exp1power';
    case 4
        FitModel = '2exp';
    case 5
        FitModel = '2exp1power';
    case 6
        FitModel = '1exp1second';
    case 7
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
    case '1exp'         % (1 exponential parameters)
        numParam = 4;
    case '1power'       % (1 exp, 1 power law parameters)
        numParam = 5;
    case '1exp1power'   % (1 exp, 1 power law parameters)
        numParam = 7;
    case '2exp'         % (2 exponential parameters)
        numParam = 6;
    case '2exp1power'   % (2 exp, 1 power law parameters)
        numParam = 8;
    case '1exp1second'  % (1 exp, 1 second order components)
        numParam = 6;
    case '2exp1second'  % (2 exp, 1 second order components)
        numParam = 8;
end

%% Restrict to fit window

TimeZero = Time(max(Data)==Data);
ConvPad = 15;    %Editable. Define how much time to add before and after the window to get rid of the edge effects of the convolution.

WindowData = Data(Time>WindowLB-ConvPad & Time<WindowUB+ConvPad);
WindowIRF = IRF(Time>WindowLB-ConvPad & Time<WindowUB+ConvPad);
WindowTime = Time(Time>WindowLB-ConvPad & Time<WindowUB+ConvPad);

%% Setup initial parameters. 
%Editable for initial guesses and bounds.

switch FitModel
    case '1exp'         % (1 exponential parameters)
        %param(1);      %Y offset
        %param(2);      %Amplitude of the exponential component
        %param(3);      %Exponential component lifetime
        %param(4);      %IRF time shift

        InitialGuess = [mean(Data(WindowBaseline)); max(Data); 5; 0];
        LowerBound = [0; 0; 0; -3];
        UpperBound = [mean(Data(WindowBaseline))*5; inf; inf; 3];

        TRPLfitfun = @(param)Fit_IRF_1exp(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

    case '1power'       % (1 power law parameters)

        %param(1);      %Y offset
        %param(2);      %Amplitude of the Power law component
        %param(3);      %Power law component onset parameter
        %param(4);      %Power law exponent parameter (alpha)
        %param(5);      %IRF time shift

        InitialGuess = [mean(Data(WindowBaseline)); max(Data)/2; 0.1; 0.5; 0];
        LowerBound = [0; 0; 0; 0; -3];
        UpperBound = [mean(Data(WindowBaseline))*5; inf; inf; 5; 3];

        TRPLfitfun = @(param)Fit_IRF_1power(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

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

        TRPLfitfun = @(param)Fit_IRF_1exp1power(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

    case '2exp'         % (2 exponential parameters)

        %param(1);      %Y offset
        %param(2);      %Amplitude of the first exponential component
        %param(3);      %First exponential component lifetime
        %param(4);      %Amplitude of the second exponential component
        %param(5);      %Second exponential component lifetime
        %param(6);      %IRF time shift

        InitialGuess = [mean(Data(WindowBaseline)); max(Data)/2; 0.5 ; max(Data)/2; 50; 0];
        LowerBound = [0; 0; 0; 0; 0; -3];
        UpperBound = [mean(Data(WindowBaseline))*5; inf; max(Time)/2; inf; 2*max(Time); 3];

        TRPLfitfun = @(param)Fit_IRF_2exp(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

    case '2exp1power'   % (2 exp, 1 power law parameters)

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
        UpperBound = [mean(Data(WindowBaseline))*5; inf; max(Time)/2; inf; 2*max(Time); inf; 0.89; 1.76; 3];

        TRPLfitfun = @(param)Fit_IRF_2exp1power(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

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

        TRPLfitfun = @(param)Fit_IRF_1exp1second(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);

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

        TRPLfitfun = @(param)Fit_IRF_2exp1second(param,WindowTime,WindowData,WindowIRF,TimeZero,ConvPad);
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
ms.XTolerance = 1E-8;   
ms.FunctionTolerance = 1E-8;
numMultiStartRuns = 50;

%gs = GlobalSearch(ms);
%rng(14,'twister')   % for reproducibility
tic                 % for timing
%[paramgs,fvalgs] = run(gs,problem);
[paramgs,fvalgs] = run(ms,problem,numMultiStartRuns);
toc


%% Generate best fit data

switch FitModel
    case '1exp'
        [FinalFit,FitYdataExp] = Decay_IRF_1exp(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
    case '1power'
        [FinalFit, FitYdataPower] = Decay_IRF_1power(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
    case '1exp1power'
        [FinalFit, FitYdataExp, FitYdataPower] = Decay_IRF_1exp1power(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
    case '2exp'
        [FinalFit, FitYdataExp, FitYdataExp2] = Decay_IRF_2exp(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
    case '2exp1power'
        [FinalFit, FitYdataExp, FitYdataExp2, FitYdataPower] = Decay_IRF_2exp1power(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
    case '1exp1second'
        [FinalFit, FitYdataExp, FitYdataSecond] = Decay_IRF_1exp1second(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
    case '2exp1second'
        [FinalFit, FitYdataExp, FitYdataExp2, FitYdataSecond] = Decay_IRF_2exp1second(paramgs,WindowTime,WindowIRF,TimeZero,ConvPad);
end

% Calculate the Chi2 from within the fit window, which is what is outputted
% by the Decay functions.Keep in mind that the WindowTime, WindowIRF, and
% WindowData have 'padding' defined by ConvPad.

KeepInd = FinalFit~=0;                                              % Find nonzero elements as can't divide by 0 to calculate the Chi-squared.
CleanTime = WindowTime(WindowTime > WindowTime(1)+ConvPad & WindowTime < WindowTime(end)-ConvPad);
CleanTime = CleanTime(KeepInd);
CleanData = WindowData(WindowTime > WindowTime(1)+ConvPad & WindowTime < WindowTime(end)-ConvPad);
CleanData = CleanData(KeepInd);
CleanFinalFit = FinalFit(KeepInd);

FitDelta = CleanData - CleanFinalFit;
FitRes = FitDelta./sqrt(CleanData);
ChiSq = sum(FitDelta.^2./CleanFinalFit);
RedChiSq = ChiSq/(length(CleanData)-length(InitialGuess)-1);        % Length Initial Guess gives you the number of fitting parameters.


%% Plot ouputs
subplot(10,1,1:6)

semilogy(Time,Data)
hold on
semilogy(CleanTime,CleanFinalFit,'LineWidth',0.75)
yline(paramgs(1),'--k')
ylabel('Counts')
fontsize(gca,12,"points")
switch FitModel
    case '1exp'
        plot(CleanTime,FitYdataExp)
        legend({'Data','Total Fit','Baseline','Exponential comp.'},...
            'Box','off','FontSize',8,'NumColumns',2)

    case '1power'
        plot(CleanTime,FitYdataPower)
        legend({'Data','Total Fit','Baseline','Power law comp.'},...
            'Box','off','FontSize',8,'NumColumns',2)

    case '1exp1power'
        plot(CleanTime,FitYdataExp)
        plot(CleanTime,FitYdataPower)
        legend({'Data','Total Fit','Baseline','Exponential comp.',...
            'Power law comp.'},'Box','off','FontSize',8,'NumColumns',2)

    case '2exp'
        plot(CleanTime,FitYdataExp)
        plot(CleanTime,FitYdataExp2)

        legend({'Data','Total Fit','Baseline','Exponential comp. 1',...
            'Exponential comp. 2'},'Box','off','FontSize',8,'NumColumns',2)
    case '2exp1power'
        plot(CleanTime,FitYdataExp)
        plot(CleanTime,FitYdataExp2)
        plot(CleanTime,FitYdataPower)
        legend({'Data','Total Fit','Baseline','Exponential comp. 1',...
            'Exponential comp. 2','Power law comp.'},'Box','off','FontSize',8,'NumColumns',2)

    case '1exp1second'
        plot(CleanTime,FitYdataExp)
        plot(CleanTime,FitYdataSecond)
        legend({'Data','Total Fit','Baseline','Exponential comp.',...
            'Second ord. comp.'},'Box','off','FontSize',8,'NumColumns',2)

    case '2exp1second'
        plot(CleanTime,FitYdataExp)
        semilogy(CleanTime,FitYdataExp2)
        semilogy(CleanTime,FitYdataSecond)
        legend({'Data','Total Fit','Baseline','Exponential comp. 1',...
            'Exponential comp. 2','Second ord. comp.'},'Box','off','FontSize',8,'NumColumns',2)

end
xlim([WindowLB WindowUB])

hold off
subplot(10,1,8:10)
plot(CleanTime,FitRes)
yline(0)
xlim([WindowLB WindowUB])
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
SearchTime = CleanTime(CleanTime>TimeZero);
SearchData = CleanFinalFit(CleanTime>TimeZero);

HalfIntTimes = SearchTime(min(abs(SearchData-max(SearchData)/2))==abs(SearchData-max(SearchData)/2));
RoughHalfIntTime = mean(HalfIntTimes);              %In case there are multiple values

DataSpacing = CleanTime(end)-CleanTime(end-1);
FineTime = RoughHalfIntTime-2*DataSpacing:DataSpacing/50:RoughHalfIntTime+2*DataSpacing;
FineData = interp1(CleanTime,CleanTime,FineTime);
FinalHalfIntTime = FineTime(min(abs(FineData-max(SearchData)/2))==abs(FineData-max(SearchData)/2));

t50 = FinalHalfIntTime(1) - TimeZero;              %TimeZero is defined as the time where the Ydata is max


% Find relative emission from the different components
switch FitModel
    case '1exp'
        RelEmExp = 1;           % Trivial since there is only 1 component.

    case '1power'
        RelEmPower = 1;         % Trivial since there is only 1 component.

    case '1exp1power' 
        ExpCompIntegral = trapz(CleanTime,FitYdataExp);
        PowerCompIntegral = trapz(CleanTime,FitYdataPower);

        RelEmExp = ExpCompIntegral/(ExpCompIntegral+PowerCompIntegral);
        RelEmPower = PowerCompIntegral/(ExpCompIntegral+PowerCompIntegral);

    case '2exp'
        Exp1CompIntegral = trapz(CleanTime,FitYdataExp);
        Exp2CompIntegral = trapz(CleanTime,FitYdataExp2);

        RelEmExp1 = Exp1CompIntegral/(Exp1CompIntegral+Exp2CompIntegral);
        RelEmExp2 = Exp2CompIntegral/(Exp1CompIntegral+Exp2CompIntegral);

    case '2exp1power'
        Exp1CompIntegral = trapz(CleanTime,FitYdataExp);
        Exp2CompIntegral = trapz(CleanTime,FitYdataExp2);
        PowerCompIntegral = trapz(CleanTime,FitYdataPower);

        RelEmExp1 = Exp1CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);
        RelEmExp2 = Exp2CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);
        RelEmPower = PowerCompIntegral/(Exp1CompIntegral+Exp2CompIntegral+PowerCompIntegral);

    case '1exp1second'
        ExpCompIntegral = trapz(CleanTime,FitYdataExp);
        SecondCompIntegral = trapz(CleanTime,FitYdataSecond);

        RelEmExp = ExpCompIntegral/(ExpCompIntegral+SecondCompIntegral);
        RelEmSecond = SecondCompIntegral/(ExpCompIntegral+SecondCompIntegral);

    case '2exp1second'

        Exp1CompIntegral = trapz(CleanTime,FitYdataExp);
        Exp2CompIntegral = trapz(CleanTime,FitYdataExp2);
        SecondCompIntegral = trapz(CleanTime,FitYdataSecond);

        RelEmExp1 = Exp1CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);
        RelEmExp2 = Exp2CompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);
        RelEmSecond = SecondCompIntegral/(Exp1CompIntegral+Exp2CompIntegral+SecondCompIntegral);

end

%% Setup fit function (1 exponential)

function RedChiSq = Fit_IRF_1exp(param,Time,Data,IRF,TimeZero,ConvPad)
% monoexponential + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
IRFshift = param(4);        %Time shift of IRF


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

% Calculate output

GuessKin = Y0 + YdataExp;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (1 power law)

function RedChiSq = Fit_IRF_1power(param,Time,Data,IRF,TimeZero,ConvPad)
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

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (1 exp, 1 power law)

function RedChiSq = Fit_IRF_1exp1power(param,Time,Data,IRF,TimeZero,ConvPad)
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

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (2 exponentials)

function RedChiSq = Fit_IRF_2exp(param,Time,Data,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
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

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
YdataExp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Calculate output

GuessKin = Y0 + YdataExp + YdataExp2;

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (2 exp, 1 power law)

function RedChiSq = Fit_IRF_2exp1power(param,Time,Data,IRF,TimeZero,ConvPad)
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

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (1 exp, 1 second order)

function RedChiSq = Fit_IRF_1exp1second(param,Time,Data,IRF,TimeZero,ConvPad)
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

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup fit function (2 exp, 1 second order)

function RedChiSq = Fit_IRF_2exp1second(param,Time,Data,IRF,TimeZero,ConvPad)
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

% Calculate Chi squared taking into consideration the padding to get rid of
% the edge effects of the convolution.

GuessKinWindow = GuessKin(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
DataWindow = Data(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

N = length(DataWindow);
ChiSq = ((DataWindow-GuessKinWindow).^2)./abs(GuessKinWindow);
RedChiSq = sum(ChiSq)/(N-length(param)-1);

end

%% Setup decay function for getting the convolved fitted data (1 exponential)

function [TotalOutput,ExpComp] = Decay_IRF_1exp(param,Time,IRF,TimeZero,ConvPad)
% Monoexponential fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp = param(2);          %Amplitude of the exponential component
ExpTau1 = param(3);         %Exponential component lifetime
IRFshift = param(4);        %Time shift of IRF

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

% Calculate output

TotalOutput = Y0 + ExpComp;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (1 power)

function [TotalOutput,PowerComp] = Decay_IRF_1power(param,Time,IRF,TimeZero,ConvPad)
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

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
PowerComp = PowerComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (1 exp, 1 power)

function [TotalOutput,ExpComp,PowerComp] = Decay_IRF_1exp1power(param,Time,IRF,TimeZero,ConvPad)
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

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
PowerComp = PowerComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (2 exponentials)

function [TotalOutput,ExpComp,ExpComp2] = Decay_IRF_2exp(param,Time,IRF,TimeZero,ConvPad)
% Two exponentials + power law fit with experimental IRF convolution

% Set up parameters.
t = Time;
Y0 = param(1);              %Y offset
ExpAmp1 = param(2);         %Amplitude of the first exponential component
ExpTau1 = param(3);         %First exponential component lifetime
ExpAmp2 = param(4);         %Amplitude of the second exponential component
ExpTau2 = param(5);         %Second exponential component lifetime
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
ExpComp = ConvYExp(1,NnegX+1:length(t)+NnegX);

% Generate Data for the second exponential component
TempYExp2 = ExpAmp2.*exp(-(t-TimeZero)./ExpTau2);
TempYExp2(t<TimeZero) = 0;
ConvYExp2 = conv(TempYExp2,ShiftedIRF/sum(ShiftedIRF));
ExpComp2 = ConvYExp2(1,NnegX+1:length(t)+NnegX);

% Calculate output

TotalOutput = Y0 + ExpComp + ExpComp2;

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp2 = ExpComp2(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (2 exp, 1 power)

function [TotalOutput,ExpComp,ExpComp2,PowerComp] = Decay_IRF_2exp1power(param,Time,IRF,TimeZero,ConvPad)
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

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp2 = ExpComp2(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
PowerComp = PowerComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (1 exp, 1 second)

function [TotalOutput,ExpComp,SecondComp] = Decay_IRF_1exp1second(param,Time,IRF,TimeZero,ConvPad)
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

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
SecondComp = SecondComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end

%% Setup decay function for getting the convolved fitted data (2 exp, 1 second)

function [TotalOutput,ExpComp,ExpComp2,SecondComp] = Decay_IRF_2exp1second(param,Time,IRF,TimeZero,ConvPad)
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

% Get output within the fit window, wihtout the padding to reduce the edge
% effects of the convolution
TotalOutput = TotalOutput(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp = ExpComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
ExpComp2 = ExpComp2(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);
SecondComp = SecondComp(Time > Time(1)+ConvPad & Time < Time(end)-ConvPad);

end