function Region = collapseLibrary(L,g,p)
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
%

% INITIALIZE
Region = struct('NBins',struct('xss',struct('xs',[]),...
                      'spectralFlux',[]),...
                 'Cell',struct('NBins',struct('xss',struct('xs',[]),...
                                     'spectralFlux',[])));

% zero-out the correct fields
for regidx = 1:g.nRegions
    for bin = [p.nGroups 110] % only the two relevant groups
        Region(regidx).NBins(L.Bin(bin)).spectralFlux = zeros(bin,1);
        for mt = L.MTs
            % region level
            Region(regidx).NBins(L.Bin(bin)).xss(L.MT(mt)).xs = [];
            
        end
        
        % cell level
        for cellidx = 1:g.regionDef(regidx).nCells
            for mt = L.MTs
                Region(regidx).Cell(cellidx).NBins(L.Bin(bin)).xss(L.MT(mt)).xs = [];
            end
            Region(regidx).Cell(cellidx).NBins(L.Bin(bin)).spectralFlux = zeros(bin,1);
        end
        
    end
end


% ACTUALLY ASSIGN STUFF NOW
for regidx = 1:g.nRegions
    R = Region(regidx);
    for cellidx = 1:g.regionDef(regidx).nCells
        RC = R.Cell(cellidx);
        cellZAIDs = g.regionDef(regidx).cellDef(cellidx).ZAIDs;
        cellAtomDen = g.regionDef(regidx).cellDef(cellidx).atomDensities; % bateman?
        cellFlux = g.regionDef(regidx).cellDef(cellidx).spectralFlux;

        % COLLAPSE CELL        
        for mt = L.MTs
            switch mt
                case 2 % scattering kernel
                    RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs = zeros(110);
                    for zaididx = 1:length(cellZAIDs)
                        RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs = RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs ...
                            + cellAtomDen(zaididx)*L.xss(L.ZAID(cellZAIDs(zaididx)),L.MT(mt)).NBins(L.Bin(110)).xs;
                    end
                case 18 % (n,f)
                    RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs = zeros(110,1);
                    for zaididx = 1:length(cellZAIDs)
                        if cellZAIDs(zaididx) > 90000 % fission only for > 90000
                            RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs = RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs ...
                                + cellAtomDen(zaididx)*L.xss(L.ZAID(cellZAIDs(zaididx)),L.MT(mt)).NBins(L.Bin(110)).xs;
                        end
                    end
                case 16 % (n,2n)
                    RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs = zeros(110,1);
                    for zaididx = 1:length(cellZAIDs)
                        if 0 %n2n exists
                            RC.xss(L.MT(mt)).NBins(L.Bin(110)).xs = RC.xss(L.MT(mt)).NBins(L.Bin(110)).xs ...
                                + cellAtomDen(zaididx)*L.xss(L.ZAID(cellZAIDs(zaididx)),L.MT(mt)).NBins(L.Bin(110)).xs;
                        end
                    end
                otherwise % all other cross sections
                    RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs = zeros(110,1);
                    for zaididx = 1:length(cellZAIDs)
                        RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs = RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs ...
                            + cellAtomDen(zaididx)*L.xss(L.ZAID(cellZAIDs(zaididx)),L.MT(mt)).NBins(L.Bin(110)).xs;
                    end                    
            end
            
            % COLLAPSE ENERGY. calls collapseEnergy for each cell
            [RC.NBins(L.Bin(p.nGroups)).xss(L.MT(mt)).xs RC.NBins(L.Bin(p.nGroups)).spectralFlux] = ...
                collapseEnergy(g.groupsIn,g.groupsOut,...
                    RC.NBins(L.Bin(110)).xss(L.MT(mt)).xs, cellFlux);
        end

        R.Cell(cellidx) = RC; %store back
    end
    
    % COLLAPSE REGION

    for mt = L.MTs
        % this is not so flexible right now. calls collapseRegion
        if g.regionDef(regidx).nCells == 2
        [R.NBins(L.Bin(p.nGroups)).xss(L.MT(mt)).xs, R.NBins(L.Bin(p.nGroups)).spectralFlux] = ...
            collapseRegion(g.regionDef(regidx).relVolumes,...
                R.Cell(1).NBins(L.Bin(p.nGroups)).xss(L.MT(mt)).xs, R.Cell(1).NBins(L.Bin(p.nGroups)).spectralFlux,...
                R.Cell(2).NBins(L.Bin(p.nGroups)).xss(L.MT(mt)).xs, R.Cell(2).NBins(L.Bin(p.nGroups)).spectralFlux);
        elseif g.regionDef(regidx).nCells == 1
        [R.NBins(L.Bin(p.nGroups)).xss(L.MT(mt)).xs, R.NBins(L.Bin(p.nGroups)).spectralFlux] = ...
            collapseRegion(g.regionDef(regidx).relVolumes,...
                R.Cell(1).NBins(L.Bin(p.nGroups)).xss(L.MT(mt)).xs, R.Cell(1).NBins(L.Bin(p.nGroups)).spectralFlux);
        end
    end

    Region(regidx) = R; %store back
    % it's weird to get a flux out from collapseRegion
end
