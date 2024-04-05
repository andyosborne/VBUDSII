function [energy mubar] = parseboundmubar(fname)
fid = fopen(fname);
if fid == -1
    error('File error');
end
energy = [];
mubar = [];
while ~feof(fid)
    fline = fgetl(fid);
    if ~isempty(strfind(fline, 'enow'))
        % Take out all alphabetic characters except E.
        fline = regexprep(fline, ' [a-zA-Z]* ', ' ');
        numbers = str2num(fline);
        energy = [energy numbers(1)];
        mubar = [mubar numbers(3)];
    end
end
end
