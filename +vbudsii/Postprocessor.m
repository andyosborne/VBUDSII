function Postprocessor(p, g, Region)

import vbudsii.*

util.PrintEntering(p, 'Postprocessor');

maxcells = 0;
for regidx = 1:g.nRegions
    if g.regionDef(regidx).nCells > maxcells
        maxcells = g.regionDef(regidx).nCells;
    end
end

if p.doPlot
    for regidx = 1:g.nRegions
        plotxsn(1,2,g.regionDef(regidx).name,{'H2O','UO2'},['k_{\infty} = ' ...
            num2str(Region(regidx).kInf)],...
            p.fineGroupDef,...
            Region(regidx).Cell(1).spectralFlux,...
            Region(regidx).Cell(2).spectralFlux);
    end
    
    figure;
    for regidx = 1:g.nRegions
        for cellidx = 1:g.regionDef(regidx).nCells
            subplot(g.nRegions,maxcells,(regidx-1)*maxcells+cellidx);
            loglog(p.fineGroupDef(1:end-1)', ...
                Region(regidx).Cell(cellidx).spectralFlux);
            title(g.regionDef(regidx).cellDef(cellidx).name);
        end
    end
end

util.PrintExiting(p, 'Postprocessor');

end
