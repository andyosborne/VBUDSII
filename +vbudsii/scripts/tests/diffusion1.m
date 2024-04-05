% 111106 this is as a result of a request to compare vbudsii, vbudsi, mcnpx and
% mcnp5.

mycd = cd;
cd(fullfile('..','..',''));
vbudsiiDir = cd;
cd(mycd); % alternatively, mfilename

%% DEFINE PARAMETER STRUCTURE
fineGroupDef = 10.^(-4:.1:7);
fewGroupDef = [1e-4 1 100e3 1e7];
p = struct('mode','diffusion'...
           'nFineGroups',length(fineGroupDef)-1,... % leave blank.
           'nFewGroups',length(fewGroupDef)-1,...
           'fineGroupDef',fineGroupDef,...
           'fewGroupDef',fewGroupDef,...
           'powerDensity',50,... % MW/m^3
           'nTimeSteps',1,...
           'vbudsiiDir',vbudsiiDir,...
           'makeLibrary',1,...
           'makeLibraryTempFlag',2,...
           'XSLibraryMAT',fullfile('..','..','data','XSLibrary.mat'),...
           'verbose',1,...
           'resolveXS',1,...
           'immutableMyMTs',1,...
           'makeRealNuFission',1,...
           'S0iterthresh',.00001,...
           'doPlot',0,...
           'diffusion',struct('nDim',1,...
                              'isBare',1,...
                              'gridpoints',300,...
                              'printDiffParams',0);

uc.pinDiam = 1; % [1,1];
uc.pinPitch = 2; % [2,2];
uc.f = 1; % some weighting of the fuel regions
uc.g = 1; % some weighting of the moderator regions
uc.sauerConst.mod = 2.35;
uc.sauerConst.fuel = 5.00;

p.uc = uc;

%{
        p = struct('nDim',1,...
                'nGroups',1,...
             'densitySet',2,...
                      'T',600,...
                     'BW',2.0246,...
             'enrichment',[0.032 0.032 0.032],...
                 'isBare',1,...
            'makeLibrary',0,...
         'printDiffParam',0,...
                     'gp',300,...
                    'BCs',[],...
              'makePlots',1,...
             'plotChoice',[1 1 1],...
               'plotSave',[0 0 0],...
               'plotName','test');
%}

%% DEFINE GEOMETRY STRUCTURE
% initialize
g = struct('nRegions',          [],...
           'regionDef',         struct('name','',...
                                       'nCells',    [],...
                                       'relVolumes',[],...
                                       'isFissile',[],...
                                       'cellDef',struct('name','',...
                                                'isFissile',[],...
                                                'initZAIDs',[],...
                                                'initNumDensities', [],...
                                                'initSpectralFlux',[]...
                                                        )));

enrichment = .032;
density_H2O = 0.72; % g/cm^3
density_UO2 = 11; %g/cm^3
[ao, wo, N_H2O] = matl([1001 2;
                        8016 1], 1, density_H2O);
[ao, wo, N_UO2] = matl([92235 enrichment;
                        92238 1-enrichment;
                        8016 2], 1, density_UO2);

% define
% one region with two cells: UO2 and H2O
g.nRegions = 4;
for regidx = 1:g.nRegions - 1
    g.regionDef(regidx).name = ['campaign' num2str(regidx)];
    g.regionDef(regidx).nCells = 2;
    g.regionDef(regidx).relVolumes = ...
        [uc.pinPitch^2-pi/4*uc.pinDiam^2, pi/4*uc.pinDiam^2]';
    % UO2 cell.
    g.regionDef(regidx).cellDef(1).name = 'UO2';
    g.regionDef(regidx).cellDef(1).initZAIDs = [92235 92238 8016];
    g.regionDef(regidx).cellDef(1).initNumDensities = N_UO2';
    g.regionDef(regidx).cellDef(1).initSpectralFlux = ones(1,p.nFineGroups);
    g.regionDef(regidx).cellDef(1).initTemp = 900;
    % H2O cell.
    g.regionDef(regidx).cellDef(2).name = 'H2O';
    g.regionDef(regidx).cellDef(2).initZAIDs = [222];
    g.regionDef(regidx).cellDef(2).initNumDensities = [N_H2O(2)];
    g.regionDef(regidx).cellDef(2).initSpectralFlux = ones(1,p.nFineGroups);
    g.regionDef(regidx).cellDef(2).initTemp = 300;
end

% Reflector region.
regrefl = 4;
g.regionDef(regrefl).name = 'reflector';
g.regionDef(regrefl).nCells = 1;
g.regionDef(regrefl).relVolumes = [1];
g.regionDef(regrefl).cellDef(1).name = 'H2O';
g.regionDef(regrefl).cellDef(1).initZAIDs = [222];
g.regionDef(regrefl).cellDef(1).initNumDensities = [N_H2O(2)];
g.regionDef(regrefl).cellDef(1).initSpectralFlux = ones(1,p.nFineGroups);
g.regionDef(regrefl).cellDef(1).initTemp = 300;

% WE STILL HAVE NO ACTUAL INPUT FOR MACROSCOPIC GEOMETRY: RADIUS OF THE
% REACTOR, SHAPE OF THE REACTOR.
% MAY WANT TO SPECIFY A REGION THAT CANNOT BE ALTERED BY BATEMAN, SUCH AS WATER
% REFLECTOR.

%{
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
%}

%% RUN VBUDSII
%cd(fullfile('..','..',''));
%addpath(fullfile('..','..',''));
addpath(p.vbudsiiDir);

tallyxs = TallyPull2('mcnpx/xs1/outp');
%p.useMCNPXTallyXS = tallyxs;

[Results, p , g, L] = main(p,g);

r = Results.Region(1);

vbudsiiPower = sum(200* r.Cell(2).fine(L.MT(18)).value .* ...
r.Cell(2).spectralFlux)

%% CHECK RESULTS
doPlot = 0;

Tally1 = TallyPull2('mcnp5/2/outp');
mcnpflux = [Tally1{2}.value{1}(1:end-1) Tally1{1}.value{1}(1:end-1)];
mcnp5flux = mcnpflux*max(max([r.Cell(1).spectralFlux r.Cell(2).spectralFlux]))/max(max(mcnpflux));

mcnp5Power = Tally1{3}.value{1}(end)

Tally1 = TallyPull2('mcnpx/2/outp');
mcnpflux = [Tally1{2}.value{1}(1:end-1) Tally1{1}.value{1}(1:end-1)];
mcnpxflux = mcnpflux*max(max([r.Cell(1).spectralFlux r.Cell(2).spectralFlux]))/max(max(mcnpflux));

mcnpxPower = Tally1{3}.value{1}(end)

load benchmarkvbuds1.mat
vbudsflux = VBUDS_flux(end:-2:2,[1 2]);
vbudsflux = vbudsflux*max(max([r.Cell(1).spectralFlux r.Cell(2).spectralFlux]))/max(max(vbudsflux));


for regidx = 1:g.nRegions
    for cellidx = 1:g.regionDef(regidx).nCells
    plotxsn(0,2,...
    [g.regionDef(regidx).name ': '...
    g.regionDef(regidx).cellDef(cellidx).name],...
    {['vbudsii. k_{\infty} = ' num2str(r.kInf)],...
     ['vbudsi, k_{\infty} = '  num2str(six_factors(7))],...
        'mcnp5. k_{eff} = 1.44','mcnpx. k_{eff} = 1.44'},...
    '',...
        p.fineGroupDef,...
        Results.Region(regidx).Cell(cellidx).spectralFlux,...
        vbudsflux(:,cellidx),...
        mcnp5flux(:,cellidx),...
        mcnpxflux(:,cellidx));
    end
end

%{
if sum(r.Cell(1).spectralFlux ~= mcnpflux(:,1)) || ...
    sum(r.Cell(2).spectralFlux ~= mcnpflux(:,2))
    disp('benchmarkmcnpx flux test FAILED');
    figure;
    subplot(1,2,1);
    loglog(p.fineGroupDef(1:end-1)'*[1 1],[r.Cell(1).spectralFlux ...
        mcnpflux(:,1)])
    title('H2O')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('vbudsii','mcnpx','Location','Best')
    subplot(1,2,2);
    loglog(p.fineGroupDef(1:end-1)'*[1 1],[r.Cell(2).spectralFlux ...
        mcnpflux(:,2)])
    title('UO2')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('vbudsii','mcnpx','Location','Best')
%   [r.Cell(2).spectralFlux./mcnpflux(:,2) r.Cell(2).spectralFlux-mcnpflux(:,2)]
    %print('benchmarkmcnp5_111026.eps','-depsc','-r300');
else
    disp('benchmarkmcnpx flux test passed');
%    [Region(1).Cell(1).spectralFlux == phi(:,1), ...
%    Region(1).Cell(2).spectralFlux == phi(:,2)]
end
%}
%{
plotxsn(1,2,g.regionDef(1).name,...
    {'vbudsi H2O','vbudsi UO2','vbudsii H2O','vbudsii UO2'},'',...
    p.fineGroupDef,...
    vbudsflux(:,2),...
    r.Cell(1).spectralFlux,...
    vbudsflux(:,1),...
    r.Cell(2).spectralFlux);


plotxsn(1,2,g.regionDef(1).name,...
    {'vbudsi H2O','vbudsi UO2','vbudsii H2O','vbudsii UO2'},'',...
    p.fineGroupDef,...
    VBUDS_flux3(:,2),...
    r.Cell(1).spectralFlux,...
    VBUDS_flux3(:,1),...
    r.Cell(2).spectralFlux);
%}
