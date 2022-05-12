%Import the raw A and B trace data for CH0 from a TDMS file.

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


%% Find and import the Raw data for CH0.

CH0GroupIndex = find(strcmp(ImportTDMS.groupNames,'CH0'));

RawAChanIndex = find(strcmp(ImportTDMS.chanNames{1,CH0GroupIndex},'A-Test  raw'));
RawBChanIndex = find(strcmp(ImportTDMS.chanNames{1,CH0GroupIndex},'B-Test  raw'));

RawADataIndex = ImportTDMS.chanIndices{1,CH0GroupIndex}(RawAChanIndex);
RawBDataIndex = ImportTDMS.chanIndices{1,CH0GroupIndex}(RawBChanIndex);

RawA = transpose(ImportTDMS.data{1,RawADataIndex});
RawB = transpose(ImportTDMS.data{1,RawBDataIndex});

%% Generate the index for the background level

BkgGroupIndex = find(strcmp(ImportTDMS.groupNames,'Background level'));
AvgBkgChanIndex = find(strcmp(ImportTDMS.chanNames{1,BkgGroupIndex},'Avg bkg level (V)'));
AvgBkgDataIndex = ImportTDMS.chanIndices{1,BkgGroupIndex}(AvgBkgChanIndex);

IndBkgChanIndex = find(strcmp(ImportTDMS.chanNames{1,BkgGroupIndex},'Individual bkg levels (V)'));
IndBkgDataIndex = ImportTDMS.chanIndices{1,BkgGroupIndex}(IndBkgChanIndex);

disp(['The background level is ' num2str(ImportTDMS.data{1,AvgBkgDataIndex}) ' V.'])
