%Import the CH0 %Abs directly from a TDMS file.

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

%% Find and import the final A-B abs data for CH0.

CH0GroupIndex = find(strcmp(ImportTDMS.groupNames,'CH0'));
AbsChanIndex = find(strcmp(ImportTDMS.chanNames{1,CH0GroupIndex},'final A-B'));
AbsDataIndex = ImportTDMS.chanIndices{1,CH0GroupIndex}(AbsChanIndex);

Abs = transpose(ImportTDMS.data{1,AbsDataIndex});