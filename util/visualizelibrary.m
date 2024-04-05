function visualizelibrary(L)
colordef white;
for z = L.ZAIDs
    for m = L.MTs
        if plotMT(z, m)
            axh = plotXStableMATLAB(m, L.groupDef(2:end), ...
                L.Ts, L.z(L.ZAID(z)).m(L.MT(m)).xs);
            title(axh, sprintf('zaid %i mt %i', z, m));
        end
    end
end
end

function truefalse = plotMT(z, m)
    if m == 102 % || m == 2
        truefalse = true;
    %if z < 90 && (m == 8 || m == 452)
    %    truefalse = false;
    %elseif m == 1 || m == 8 || m == 102 || m == 251 || m == 452
    %    truefalse = true;
    else
        truefalse = false;
    end
%    if m == 2 || (m >= 50 && m <=91)
%        truefalse = false;
%    else
%        truefalse = true;
%    end
end

function axh = plotXStableMATLAB(m, xvals, Ts, table)

col0rs = [79 129 189;
          192 80 77;
          155 187 89;
          128 100 162;
          247 150 70;
          0 176 240;
          0 0 255;
         0 128 0;
         255 0 0;
         0 191 191;
         191 0 191;
         191 191 0;
         64 64 64];


         [imin imax] = energy2idxrange(xvals, -1, 10e10);

    figh= figure;
    axh = axes;
    hold on;
    templegend = cell(1, length(Ts));
    for iT = 1:length(Ts)
        if m == 2
        fillh = fill(log([xvals(imin:imax) xvals(imax:-1:imin)]), ...
            log([sum(table(:,imin:imax,iT,1))'; sum(table(:,imax:-1:imin,iT,end))']), ...
            col0rs(iT,:)/255, ...
            'EdgeColor', col0rs(iT,:)/255 ...
            );
        else
        fillh = fill(log([xvals(imin:imax) xvals(imax:-1:imin)]), ...
            log([table(imin:imax,iT,1); table(imax:-1:imin,iT,end)]), ...
            col0rs(iT,:)/255, ...
            'EdgeColor', col0rs(iT,:)/255 ...
            );
    %        'EdgeColor', 'none' ...
        end
        alpha(fillh, .1);
        templegend{iT} = num2str(Ts(iT));
    end
    legend(templegend);
    % TODO temperature.
    xlabel('energy (eV)');
    ylabel('\sigma (b)');
end


function [imin imax] = energy2idxrange(evals, emin, emax)

    imin = max(1, find( emin < evals, 1, 'first'));
    imax = min(length(evals), find( emax >= evals, 1, 'last'));
    
end
