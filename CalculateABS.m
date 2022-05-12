%Calculate %Abs from the raw data in a TDMS file.

%Reflection mode calculation.

%% Point to the file.

[FileTDMS, PathTDMS] = uigetfile('D:\TAM data\Robert\*.tdms');

ImportTDMS = TDMS_readTDMSFile([PathTDMS, FileTDMS]);

%% Find where the time data is and import it.

TimeGroupIndex = find(strcmp(ImportTDMS.groupNames,'Time'));
TimeDataIndex = ImportTDMS.chanIndices{1,TimeGroupIndex};

%Check the first 10 index differences to get a delta double array
delta = ImportTDMS.data{1,TimeDataIndex}(2:11) - ImportTDMS.data{1,TimeDataIndex}(1:10);
int_delta = int16(delta);

%Check if the double array has equal integer values to distinguish time format
if delta == int_delta
    Time = transpose(ImportTDMS.data{1,TimeDataIndex} * 4e-9);
else
    Time = transpose(ImportTDMS.data{1,TimeDataIndex});
end


%% Find the Raw data for CH0, as well as the background voltage, and calculate the abs with it.

CH0GroupIndex = find(strcmp(ImportTDMS.groupNames,'CH0'));

RawAChanIndex = find(strcmp(ImportTDMS.chanNames{1,CH0GroupIndex},'A-Test  raw'));
RawBChanIndex = find(strcmp(ImportTDMS.chanNames{1,CH0GroupIndex},'B-Test  raw'));

RawADataIndex = ImportTDMS.chanIndices{1,CH0GroupIndex}(RawAChanIndex);
RawBDataIndex = ImportTDMS.chanIndices{1,CH0GroupIndex}(RawBChanIndex);

BkgGroupIndex = find(strcmp(ImportTDMS.groupNames,'Background level'));
BkgChanIndex = find(strcmp(ImportTDMS.chanNames{1,BkgGroupIndex},'Avg bkg level (V)'));
BkgDataIndex = ImportTDMS.chanIndices{1,BkgGroupIndex}(BkgChanIndex);

Gain = 405; %Change gain as appropriate
%Abs = -(RawA - RawB)/Background. When the active filter isn't used and the system is DC coupled, use %Abs = -(RawA - RawB - Background)/Background
Abs = transpose(-(ImportTDMS.data{1,RawADataIndex}-ImportTDMS.data{1,RawBDataIndex}/(ImportTDMS.data{1,BkgDataIndex}*Gain))); %For CH0

disp(['The background level is ' num2str(ImportTDMS.data{1,BkgDataIndex}) ' V.'])
