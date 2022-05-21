% Pointing to a certain directory with a specific file type:
clear % To ensure workspace refreshes each time script is run
path = uigetdir();
FolderContent = dir(path);

file = {};
for FileIndex = 1:1:size(FolderContent,1)
    if endsWith(FolderContent(FileIndex).name, '.txt') == 1
        file{end+1,1} = readmatrix([path,'\',FolderContent(FileIndex).name]);

        hold on % This line allows to plot multiple lines on the same axis
        line = plot(file{end,1}(:,1),file{end,1}(:,2));
        line.DisplayName = FolderContent(FileIndex).name; % This assigns name of line object in legend
    end
end
hold off % For consistency turn hold off
legend("Box","off"); % You can change this to "on" for a box around legend
xlim("tight") % You can change this to "padded" for some space around data
ylim("tight") % You can change this to "padded" for some space around data