function g = makeGeometry(p)
%MAKEGEOMETRY Create a geometry structure from a parameter structure p.
% Description: The geometry structure defines the cells of regions, and the
% materials in each cell for each region. Also holds the flux for each
% cell. This function does not create a general geometry structure, but
% rather a specific one, with parameters filled in. This means that a new
% MAKEGEOMETRY function would be written if the structure g is to change.
%
% USE:  g = makeGeometry(p)
%
% NOTES: The particular geometry here is 3 regions, all the same, with two
% cells: UO2 and H2O. A fourth region is beyond the reactor core, and is
% here defined to be a water reflector region. This fourth region can be
% overwritten later on by modifying diffusion parameters.
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110519 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%
% DEPENDENCIES: 
%   FLUXDATA.mat
%   reducedlib.mat
%

% load mat files: holds flux for each cell, and raw cross sections
load FLUXDATA.mat;
load reducedlib.mat;

% initialize output
g = struct('nRegions',4,... % 4th region is reflector or vacuum. 
    'regionDef',struct('nCells',[],...
    'relVolumes',[],...
    'cellDef',struct('ZAIDs',[],...
    'atomDensities',[],...
    'spectralFlux',[])),...
    'groupsIn',edges,...
    'groupsOut',[]);

% REGION-CELL DEFINITIONS
for regidx = 1:g.nRegions-1
    g.regionDef(regidx).nCells = 2;
    g.regionDef(regidx).relVolumes = [0.5271 0.9493]';
    
    % UO2 cell
    g.regionDef(regidx).cellDef(1).ZAIDs = [92235 92238 8016];
    g.regionDef(regidx).cellDef(1).atomDensities = makeUO2(p.densitySet,p.enrichment(regidx),'macro'); % needs to know temperature
    g.regionDef(regidx).cellDef(1).spectralFlux = FLUX(219:-2:1,1);
    
    % H2O cell
    g.regionDef(regidx).cellDef(2).ZAIDs = [222]; % just water
    g.regionDef(regidx).cellDef(2).atomDensities = makeH2O(p.densitySet,'macro'); % needs to know temperature
    g.regionDef(regidx).cellDef(2).spectralFlux = FLUX(219:-2:1,2);
end

% reflector region
regrefl = 4;
g.regionDef(regrefl).nCells = 1;
g.regionDef(regrefl).relVolumes = [1];
g.regionDef(regrefl).cellDef(1).ZAIDs = [222];
g.regionDef(regrefl).cellDef(1).atomDensities = makeH2O(p.densitySet,'macro');
g.regionDef(regrefl).cellDef(1).spectralFlux = FLUX(219:-2:1,2);

if p.nGroups == 3
    g.groupsOut = [1e-4 1 100e3 10e6];
elseif p.nGroups == 1
    g.groupsOut = [1e-4 10e6];
end
end