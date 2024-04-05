function str = tallystring(varargin)

inps = varargin;
str = num2str(inps{1});
for i = 2:length(inps)
    str = [str '.' num2str(inps{i})];


end

