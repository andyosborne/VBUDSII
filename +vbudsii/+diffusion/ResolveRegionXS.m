function Region = ResolveRegionXS(L, p, g, Region)
%COLLAPSELIBRARY Collapse the cross section Library by cell definitions, by
%neutron energy, and by region definition.
% Description:  The cross sections must be prepared from the Library in a
% form suitable for diffusion. This function is a first part of that
% process. The cross section data must reflect the geometry of the reactor,
% and so the cross sections for the cells are developed, then those cross
% sections are brought down to just three energy groups using cell fluxes
% from Multi-cell, and those cross sections are smeared within each region
% to provide 3-group per-region cross sections. Accordingly, the function
% requires the geometry structure g. As well, it requires some of the
% information in the parameter structure p. Output is stored in the
% structure Region.
%
% USE:  Region = collapseLibrary(L,g,p)
%
% NOTES:
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110519 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%   1- n,2n cross section computation
%
% DEPENDENCIES:
%   collapseEnergy
%   collapseRegion
% TASKLIST:
%   1- Error-check that fewFlux is not changing for each MT.
%   2- CollapseRegion must be rewritten to manage more flexible inputs.
%

% ACTUALLY ASSIGN STUFF NOW
for regidx = 1:g.nRegions
    for cellidx = 1:g.regionDef(regidx).nCells

        for mt = L.MTs

            % COLLAPSE ENERGY.
            % Call CollapseEnergy for each cell.
            [Region(regidx).Cell(cellidx).few(L.MT(mt)).value ...
             Region(regidx).Cell(cellidx).fewFlux] = ...
                CollapseEnergy(p.fineGroupDef, p.fewGroupDef, ...
                    Region(regidx).Cell(cellidx).fine(L.MT(mt)).value, ...
                    Region(regidx).Cell(cellidx).spectralFlux);

            % Error-check that fewFlux is not changing for each MT.
            
        end

    end
    
    % COLLAPSE REGION

    for mt = L.MTs

        % This is not so flexible right now. Calls CollapseRegion.
        if g.regionDef(regidx).nCells == 2
            [Region(regidx).few(L.MT(mt)).value ...
             Region(regidx).fewFlux] = ...
                CollapseRegion(g.regionDef(regidx).relVolumes,...
                    Region(regidx).Cell(1).few(L.MT(mt)).value,...
                    Region(regidx).Cell(1).fewFlux,...
                    Region(regidx).Cell(2).few(L.MT(mt)).value,...
                    Region(regidx).Cell(2).fewFlux);
        elseif g.regionDef(regidx).nCells == 1
            [Region(regidx).few(L.MT(mt)).value ...
             Region(regidx).fewFlux] = ...
                CollapseRegion(g.regionDef(regidx).relVolumes,...
                    Region(regidx).Cell(1).few(L.MT(mt)).value,...
                    Region(regidx).Cell(1).fewFlux);
        else
            error(['ERROR: ResolveRegionXS cannot handle more than 2 '...
                'cells the way the code is currently written']);
        end
    end

end
