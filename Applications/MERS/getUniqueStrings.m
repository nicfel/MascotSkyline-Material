function [uniqueStrings, ind] = getUniqueStrings(inputCell)
% Compare the strings pairwise, ignoring the '-' character
uniqueStrings = {inputCell{1}};
ind = zeros(1,length(inputCell));
ind(1)=1;
for i = 2:numel(inputCell)
    found = false;
    for j = 1:numel(uniqueStrings)
        if length(inputCell{i})==length(uniqueStrings{j})
            diffs = find(inputCell{i} ~= uniqueStrings{j});

            gaps = sum(inputCell{i}(diffs)=='-') + sum(uniqueStrings{j}(diffs)=='-');
            if gaps==length(diffs)
                ind(i)=j;
                found = true;
                break;
            end
        end
    end
    if ~found
        uniqueStrings{end+1} = inputCell{i};
        ind(i)=length(uniqueStrings);
    end
end
end
