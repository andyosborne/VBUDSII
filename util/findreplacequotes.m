function fname_out = findreplacequotes(fname)

fname_out = [fname '_revised.m'];
fidin = fopen([fname '.m'], 'r');
fidout = fopen(fname_out, 'w');

while ~feof(fidin)
    fline = fgetl(fidin);
    quoteidxs = find(fline == '"');
    if quoteidxs
        for quoteidx = quoteidxs 
            fline(quoteidx) = '''';
        end
    end
    fprintf(fidout, [fline '\n']);
end

fclose(fidin);
fclose(fidout);

end

