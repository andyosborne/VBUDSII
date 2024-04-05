function spectralFlux = extractmcnpflux(outpdirname)

spectralFlux = zeros(110,2);
fineGroupDef = zeros(110,1);

fid = fopen(outpdirname,'r');
L = fgets(fid);
linenum = 1;
cellidx = 0;
while L~=-1 % end of file
    linenum = linenum + 1;
    if ~isempty(strfind(L, '  fluxcapacitor'))
        while L~=-1
            if ~isempty(strfind(L, 'energy'))
                L = fgets(fid);
                linenum = linenum + 1;
                cellidx = cellidx + 1
                i = 0;
                while isempty(strfind(L, 'total'));
                    i = i + 1;
                    line2str = str2num(L);
                    line2str;
                    fineGroupDef(i) = line2str(1);
                    spectralFlux(i,cellidx) = line2str(2);
                    L = fgets(fid);
                    linenum = linenum + 1;
                end
                if cellidx == 2
                    break;
                end
            end
            if cellidx == 2
                break;
            end
            L = fgets(fid);
            linenum = linenum + 1;
        end
    end
    L = fgets(fid);
    inenum = linenum + 1;
end

fclose(fid);

figure;
subplot(1,2,1);
loglog(fineGroupDef,spectralFlux(:,2))
title('H2O')
ylabel('\phi (# n cm^{-2} s^{-1})')
xlabel('E (eV)')
subplot(1,2,2);
loglog(fineGroupDef,spectralFlux(:,1))
title('UO2')
ylabel('\phi (# n cm^{-2} s^{-1})')
xlabel('E (eV)')
%   print('benchmarkvbuds1_111014.eps','-depsc','-r300');
end
