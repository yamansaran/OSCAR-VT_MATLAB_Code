fid = fopen('Planet4589Active.txt', 'r');
lines = {};
while ~feof(fid)
    line = fgetl(fid);
    if ~startsWith(strtrim(line), '#') && ~isempty(strtrim(line))
        lines{end+1} = line;
    end
end
fclose(fid);

% Split on 2+ whitespace
header = regexp(lines{1}, '\s{2,}', 'split');
header = strtrim(header);
header = header(~cellfun('isempty', header));  % remove empty cells

numCols = length(header);

data = cell(length(lines)-1, numCols);
for i = 2:length(lines)
    row = regexp(lines{i}, '\s{2,}', 'split');
    row = strtrim(row);
    row = row(~cellfun('isempty', row));  % remove empty cells
    
    % Match row length to header length
    if length(row) > numCols
        row = row(1:numCols);  % truncate if too many
    elseif length(row) < numCols
        row(end+1:numCols) = {''};  % pad if too few
    end
    
    data(i-1, :) = row;
end

% Create table and write
T = cell2table(data, 'VariableNames', matlab.lang.makeValidName(header));
writetable(T, 'Planet4589Objects.csv');