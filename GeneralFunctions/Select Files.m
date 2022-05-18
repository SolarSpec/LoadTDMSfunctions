% Select certain files to import and plot
clear % To ensure workspace refreshes each time script is run
[FolderContent,path] = uigetfile('*.txt;*.csv;*.xlsx','Select one or more input files','MultiSelect','on');
FolderContent = FolderContent';
data = {};
for FileIndex = 1:1:size(FolderContent,2)
        data{FileIndex,1} = readmatrix(string(strjoin([path,FolderContent(FileIndex)],"")));

        hold on % This line allows to plot multiple lines on the same axis
        line = plot(data{FileIndex,1}(:,1),data{FileIndex,1}(:,2));
        line.DisplayName = string(FolderContent(FileIndex)); % This assigns name of line object in legend
end
hold off % For consistency turn hold off
legend("Box","off"); % You can change this to "on" for a box around legend
xlim("tight") % You can change this to "padded" for some space around data
ylim("tight") % You can change this to "padded" for some space around data