% Find maximum value and index of matrix
challenge = magic(5); % Create 5x5 matrix with arbitrary values
MaxRows = max(challenge); % Finds maximum of each row
MaxVal = max(max(challenge)); % Finds total maximum of entire matrix

linInd = find(challenge==MaxVal); % Finds linear index of maximum element of entire matrix
[m,n] = ind2sub(linInd); % Finds row (m) and column (n) of the certain index in question

%% More efficient commands
% Find maximum value and index of same matrix using one line
[Mrows,ICol] = max(challenge); % Finds maximum values of each row while giving index of each row per column
[MRows,ILin] = max(challenge, [],'linear'); % Finds maximum values of each row while giving linear indices of each
[M,I] = max(challenge,[],'all'); % Finds total maximum of entire matrix with linear index max element
% You can again use [m,n] = ind2sub(I) to get this linear index into rows and columns
