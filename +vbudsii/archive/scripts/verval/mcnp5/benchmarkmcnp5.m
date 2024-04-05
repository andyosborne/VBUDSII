% This tests Multicell only against Geoff's original code.
% must specify exactly what this tests: tests the modified Multicell(geoff's
% vbudsii) against Geoff's original VBUDSII. Not a test for physical accuracy.

% do some cd-funk to get the absolute path of the vbudsii installation, then cd
% back into the directory that this file is in.

% MULTICELL4 (110918)

mycd = cd;
cd(fullfile('..','..','..',''));
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
           'makeLibraryTempFlag',2,...
           'XSLibraryMAT',fullfile('..','..','data','XSLibrary.mat'),...
           'verbose',1,...
           'resolveXS',1,...
           'immutableMyMTs',1,...
           'makeRealNuFission',1,...
           'S0iterthresh',.00001); % this largely affects the validation effort.

uc.pinDiam = 1; % [1,1];
uc.pinPitch = 2; % [2,2];
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
g.nRegions = 1;
g.regionDef(1).name = 'campaign1';
g.regionDef(1).nCells = 2;
g.regionDef(1).relVolumes = [uc.pinPitch^2-pi/4*uc.pinDiam^2, pi/4*uc.pinDiam^2];
% H2O cell
g.regionDef(1).cellDef(1).name = 'H2O';
g.regionDef(1).cellDef(1).initZAIDs = [222];
g.regionDef(1).cellDef(1).initNumDensities = [N_H2O(2)];
g.regionDef(1).cellDef(1).initSpectralFlux = ones(1,p.nFineGroups);
% UO2 cell
g.regionDef(1).cellDef(2).name = 'UO2';
g.regionDef(1).cellDef(2).initZAIDs = [92235 92238 8016];
g.regionDef(1).cellDef(2).initNumDensities = N_UO2';
g.regionDef(1).cellDef(2).initSpectralFlux = ones(1,p.nFineGroups);

%% RUN VBUDSII
%cd(fullfile('..','..',''));
%addpath(fullfile('..','..',''));
addpath(p.vbudsiiDir);

[Results, p , g] = main(p,g);

r = Results.Region(1);

%% CHECK RESULTS
doPlot = 0;

Tally1 = TallyPull2('2/outp');
mcnpflux = [Tally1{2}.value{1}(1:end-1) Tally1{1}.value{1}(1:end-1)];
mcnpflux = mcnpflux*max(max([r.Cell(1).spectralFlux r.Cell(2).spectralFlux]))/max(max(mcnpflux));

if sum(r.Cell(1).spectralFlux ~= mcnpflux(:,1)) || ...
    sum(r.Cell(2).spectralFlux ~= mcnpflux(:,2))
    disp('benchmarkmcnp5 flux test FAILED');
    figure;
    subplot(1,2,1);
    loglog(p.fineGroupDef(1:end-1)'*[1 1],[r.Cell(1).spectralFlux ...
        mcnpflux(:,1)])
    title('H2O')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('vbudsii','mcnp5','Location','Best')
    subplot(1,2,2);
    loglog(p.fineGroupDef(1:end-1)'*[1 1],[r.Cell(2).spectralFlux ...
        mcnpflux(:,2)])
    title('UO2')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('vbudsii','mcnp5','Location','Best')
    [r.Cell(2).spectralFlux./mcnpflux(:,2) r.Cell(2).spectralFlux-mcnpflux(:,2)]
    %print('benchmarkmcnp5_111026.eps','-depsc','-r300');
else
    disp('benchmarkmcnp5 flux test passed');
%    [Region(1).Cell(1).spectralFlux == phi(:,1), ...
%    Region(1).Cell(2).spectralFlux == phi(:,2)]
end

% calculate total power.

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
