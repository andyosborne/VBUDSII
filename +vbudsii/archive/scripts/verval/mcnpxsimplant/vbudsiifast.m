function [mcnpxflux TallyP vbudsiiflux r L p g] = vbudsiifast(runMCNPX)

    % fast reactor for paper
    if nargin == 0
        runMCNPX = 0;
    end

    mycd = cd;

cd(fullfile('..','..',''));
vbudsiiDir = cd;
cd(mycd); % alternatively, mfilename


%% PL PL PL PL PL PL PL PL
pL = struct('metamat','FarmXSMetaP1Fast',...
            'address','/home/fitze/LANL/XSFarmP1Fast',...
            'verbose',1,...
            'explore',0,...
               'mine',0,...
              'plant',0,...
          'fertilize',0,...
            'harvest',1,...
                 'e',struct('searchWeb',0),...
                 'm',struct('downloadTape20s',0),...
                 'p',struct('MFMTs',[3 1;
                                     3 4;
                                     3 18;
                                     3 102;
                                     3 251;
                                     3 452],...
                              'MFs',[6],...
                               'Ts',[296 450 600 1000],...
                              'S0s',10.^([10 5 3 2 1 0 -1]),...
                           'IGNstr','1',...
                         'groupDef',10.^(-4:.1:7),...
                           'IWTstr','3'),...
                 'f',struct('overwriteOutput',1,...
                                 'logResults',0),...
                 'h',struct('makeLibrary',1,...
                             'makeArrays',0),...
             'ZAIDs',[1001 8016 11023 92235 92238]);
%                               'Ts',[296 350 400 450 500 600 800 1000],...
                      %        'MFs',[6],...
                %                         'getInelasticParts',1,...
                     %'farmInelastic',1,...

%% DEFINE PARAMETER STRUCTURE P P P P P P P
fineGroupDef = 10.^(-4:.1:7);
fewGroupDef = [1e-4 1 100e3 1e7];
p = struct('nFineGroups',length(fineGroupDef)-1,... % leave blank.
    'nFewGroups',length(fewGroupDef)-1,...
    'fineGroupDef',fineGroupDef,...
    'fewGroupDef',fewGroupDef,...
    'powerDensity',50,... % MW/m^3
    'nTimeSteps',1,...
    'vbudsiiDir',vbudsiiDir,...
    'makeLibrary',1,...
    'makeLibraryTempFlag',pL,... %3,...
    'XSLibraryMAT','XSLibraryFast.mat',...
    'verbose',1,...
    'resolveXS',1,...
    'immutableMyMTs',1,...
    'makeRealNuFission',1,...
    'S0iterthresh',.00001,...
    'doPlot',0); % this largely affects the validation effort.

% Determining pin diameter.
uc.pinPitch = 2; % [2,2];

fuelVolFraction = 0.37;
uc.pinDiam = sqrt(uc.pinPitch^2 * fuelVolFraction * 4 / pi);

uc.f = 1; % some weighting of the fuel regions
uc.g = 1; % some weighting of the moderator regions
uc.sauerConst.mod = 2.35;
uc.sauerConst.fuel = 5.00;

p.uc = uc;

%% DEFINE GEOMETRY STRUCTURE
% initialize
g = struct('nRegions',          [],...
    'regionDef',         struct('name','',...
    'nCells',    [],...
    'relVolumes',[],...
    'isFissionable',[],...
    'cellDef',struct('name','',...
    'isFissionable',[],...
    'initZAIDs',[],...
    'initNumDensities',[],...
    'initDensity',[],...
    'initSpectralFlux',[],...
    'initTemp',[]...
    )));
%% G G G G G G G G G G G G
enrichment = .256;
density_Na = 0.882; % g/cm^3
density_UO2 = 11; %g/cm^3
[ao, wo, N_Na] = matl([11023 1], 1, density_Na);
[ao, wo, N_UO2] = matl([92235 enrichment;
    92238 1-enrichment;
    8016 2], 1, density_UO2);

% define
% one region with two cells: UO2 and H2O
g.nRegions = 1;
g.regionDef(1).name = 'campaign1';
g.regionDef(1).nCells = 2;
g.regionDef(1).relVolumes = [uc.pinPitch^2-pi/4*uc.pinDiam^2, pi/4*uc.pinDiam^2];
% H2O cell
g.regionDef(1).cellDef(1).name = 'Na';
g.regionDef(1).cellDef(1).initZAIDs = [11023];
g.regionDef(1).cellDef(1).initNumDensities = [N_Na];
g.regionDef(1).cellDef(1).initDensity = density_Na;
g.regionDef(1).cellDef(1).initSpectralFlux = ones(1,p.nFineGroups);
g.regionDef(1).cellDef(1).initTemp = 300 + 273;
% UO2 cell
g.regionDef(1).cellDef(2).name = 'UO2';
g.regionDef(1).cellDef(2).initZAIDs = [92235 92238 8016];
g.regionDef(1).cellDef(2).initNumDensities = N_UO2';
g.regionDef(1).cellDef(2).initDensity = density_UO2;
g.regionDef(1).cellDef(2).initSpectralFlux = ones(1,p.nFineGroups);
g.regionDef(1).cellDef(2).initTemp = 900 + 273;

%% RUN MCNPX

ifname = 'fast';
mcnpxexecpath = fullfile('mcnpxfast');

if runMCNPX
    % ~exist([mcnpxexecpath filesep 'outp'], 'file')
    writeMCNPXSpectralInput(p, g, mcnpxexecpath, ifname);
    cd(mcnpxexecpath);
    system(['mcnpx i= ' ifname]);
    system('rm runtpe');
end

cd(mycd);

if runMCNPX
    [Tally1 TallyP] = TallyPull2([mcnpxexecpath filesep 'outp']);
    save mcnpxfast Tally1 TallyP
else
    load mcnpxfast
end

%% RUN VBUDSII
addpath(p.vbudsiiDir);

tallyxs = Tally1;

[Results, p , g, L] = main(p,g);

r = Results.Region(1);

cd(mycd);

%% SCALE FLUX

% vbudsii
vbudsiiflux = [Results.Region(1).Cell(1).spectralFlux  ...
                Results.Region(1).Cell(2).spectralFlux];

vbudsiiPowerDensity = r.Cell(2).powerDensity;
vbudsiiFluxIntegral = ...
    sum(diff(L.groupDef').*sum([Results.Region(1).Cell(1).spectralFlux  ...
                Results.Region(1).Cell(2).spectralFlux], 2));
vbudsiiflux = vbudsiiflux/vbudsiiPowerDensity;

% mcnpx
mcnpflux = [Tally1{2}.value{1}(1:end-1) Tally1{1}.value{1}(1:end-1)];

mcnpxPower = TallyP('14.100.-6.-8');
mcnpxPowerDensity = mcnpxPower(end);
mcnpxFluxIntegral = sum(diff(L.groupDef').*sum(mcnpflux, 2));

mcnpxflux = mcnpflux*max(max(vbudsiiflux))/max(max(mcnpflux));
mcnpxflux = mcnpflux/mcnpxPowerDensity;
%mcnpxflux = mcnpflux/mcnpxFluxIntegral;

% make MCNPX 1-group cross sections.
end
