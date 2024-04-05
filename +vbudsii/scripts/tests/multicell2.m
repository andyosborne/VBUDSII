% This tests Multicell only against Geoff's original code.
% must specify exactly what this tests: tests the modified Multicell(geoff's
% vbudsii) against Geoff's original VBUDSII. Not a test for physical accuracy.

% do some cd-funk to get the absolute path of the vbudsii installation, then cd
% back into the directory that this file is in.

mycd = cd;
cd(fullfile('..','..',''));
vbudsiiDir = cd;
cd(mycd); % alternatively, mfilename

%% DEFINE PARAMETER STRUCTURE
fineGroupDef = 10.^(-4:.1:7);
fewGroupDef = [1e-4 1 100e3 1e7];
p = struct('nFineGroups',length(fineGroupDef)-1,... % leave blank.
           'nFewGroups',length(fewGroupDef)-1,...
           'fineGroupDef',fineGroupDef,...
           'fewGroupDef',fewGroupDef,...
           'temp',600,...
           'powerDensity',50,... % MW/m^3
           'nTimeSteps',1,...
           'vbudsiiDir',vbudsiiDir,...
           'makeLibrary',1,...
           'XSLibraryMAT',fullfile('..','..','data','XSLibrary.mat'),...
           'verbose',1);

%% DEFINE GEOMETRY STRUCTURE
% initialize
g = struct('nRegions',          [],...
           'regionDef',         struct('name','',...
                                       'nCells',    [],...
                                       'relVolumes',[],...
                                       'isFissile',[],...
                                       'cellDef',struct('name','',...
                                                'initZAIDs',[],...
                                                'initNumDensities', [],...
                                                'initSpectralFlux',[]...
                                                        )));

% define
% one region with two cells: UO2 and H2O
g.nRegions = 1;
g.regionDef(1).name = 'campaign1';
g.regionDef(1).nCells = 2;
g.regionDef(1).relVolumes = [4-pi/4, pi/4];
% H2O cell
g.regionDef(1).cellDef(1).name = 'H2O';
g.regionDef(1).cellDef(1).initZAIDs = [1001 8016];
g.regionDef(1).cellDef(1).initNumDensities = [.145/3, 2/3*.145];
g.regionDef(1).cellDef(1).initSpectralFlux = ones(1,p.nFineGroups);
% UO2 cell
g.regionDef(1).cellDef(2).name = 'UO2';
g.regionDef(1).cellDef(2).initZAIDs = [92235 92238 8016];
g.regionDef(1).cellDef(2).initNumDensities = 20*[.0003867, .007347, .01547];
g.regionDef(1).cellDef(2).initSpectralFlux = ones(1,p.nFineGroups);

%% RUN VBUDSII
%cd(fullfile('..','..',''));
%addpath(fullfile('..','..',''));
addpath(p.vbudsiiDir);

Region = main(p,g);

%% CHECK RESULTS
load multicell2data.mat
if sum(Region(1).Cell(1).spectralFlux ~= phi(:,1)) || ...
    sum(Region(1).Cell(2).spectralFlux ~= phi(:,2))
    disp('multicell2 test failed');
else
    disp('multicell2 test passed');
%    [Region(1).Cell(1).spectralFlux == phi(:,1), ...
%    Region(1).Cell(2).spectralFlux == phi(:,2)]
end

%load FLUXDATA.mat;
%load reducedlib.mat;
%
%% initialize output
%g = struct('nRegions',4,... % 4th region is reflector or vacuum.
%    'regionDef',struct('nCells',[],...
%    'relVolumes',[],...
%    'cellDef',struct('ZAIDs',[],...
%    'atomDensities',[],...
%    'spectralFlux',[])),...
%    'groupsIn',edges,...
%    'groupsOut',[]);
%
%% REGION-CELL DEFINITIONS
%for regidx = 1:g.nRegions-1
%    g.regionDef(regidx).nCells = 2;
%    g.regionDef(regidx).relVolumes = [0.5271 0.9493]';
%
%    % UO2 cell
%    g.regionDef(regidx).cellDef(1).ZAIDs = [92235 92238 8016];
%    g.regionDef(regidx).cellDef(1).atomDensities = makeUO2(p.densitySet,p.enrichment(regidx),'macro'); % needs to know temperature
%    g.regionDef(regidx).cellDef(1).spectralFlux = FLUX(219:-2:1,1);
%
%    % H2O cell
%    g.regionDef(regidx).cellDef(2).ZAIDs = [222]; % just water
%    g.regionDef(regidx).cellDef(2).atomDensities = makeH2O(p.densitySet,'macro'); % needs to know temperature
%    g.regionDef(regidx).cellDef(2).spectralFlux = FLUX(219:-2:1,2);
%end
%
%% reflector region
%regrefl = 4;
%g.regionDef(regrefl).nCells = 1;
%g.regionDef(regrefl).relVolumes = [1];
%g.regionDef(regrefl).cellDef(1).ZAIDs = [222];
%g.regionDef(regrefl).cellDef(1).atomDensities = makeH2O(p.densitySet,'macro');
%g.regionDef(regrefl).cellDef(1).spectralFlux = FLUX(219:-2:1,2);
%
%if p.nGroups == 3
%    g.groupsOut = [1e-4 1 100e3 10e6];
%elseif p.nGroups == 1
%    g.groupsOut = [1e-4 10e6];
%end
%end


%struct   double 1 x M, M = 110
%    fewGroupStruct  double 1 x Q, Q = 3
%        nRegions        int, 3
%                regionDef(regidx).  nCells          int, 3
%                            regidx= relVolumes      1 x nCells
%                                        1 outer cellDef(cellidx).   initnZAIDs
%                                        int
%                                                    2 middle        cellidx=
%                                                    initZAIDs           1 x M
%                                                                3 inner     1
%                                                                UO2
%                                                                initAtomDensities
%                                                                1 x initnZAIDs
%                                                                                        2
%                                                                                        H2O
%                                                                                        initSpectralFlux
%                                                                                        Q
%                                                                                        x
%                                                                                        initnZAIDs

























