%% Select certain files to import and plot. USUALLY SAME SIZE DATA
clear % To ensure workspace refreshes each time script is run

% Select the files and get number of files
[FolderContent,path] = uigetfile('*.txt;*.csv;*.xlsx','Select one or more input files','MultiSelect','on');
FolderContent = FolderContent';
L = length(FolderContent);

data = cell(L,1); % Preallocate cell array to hold each x and y column
for FileIndex = 1:L
        % Read X and Y data of each file
        data{FileIndex,1} = readmatrix(string(strjoin([path,FolderContent(FileIndex)],"")));

        hold on % This line allows to plot multiple lines on the same axis
        line = plot(data{FileIndex,1}(:,1),data{FileIndex,1}(:,2));

        % This assigns name of line object in legend
        line.DisplayName = string(FolderContent(FileIndex)); 
end
hold off % For consistency turn hold off
legend("Box","off"); % You can change this to "on" for a box around your legend
xlim("tight") % You can change this to "padded" for some space around your data
ylim("tight") % You can change this to "padded" for some space around your data