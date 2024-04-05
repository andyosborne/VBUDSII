

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
           'makeLibraryTempFlag',2,...
           'XSLibraryMAT',fullfile('..','..','data','XSLibrary.mat'),...
           'verbose',1,...
           'resolveXS',1,...
           'immutableMyMTs',1,...
           'makeRealNuFission',0,...
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
g.regionDef(1).relVolumes = [uc.pinPitch^2 - pi/4*uc.pinDiam^2, pi/4*uc.pinDiam^2];
% H2O cell
g.regionDef(1).cellDef(1).name = 'H2O';
g.regionDef(1).cellDef(1).initZAIDs = [222];
g.regionDef(1).cellDef(1).initNumDensities = [N_H2O(2)];
g.regionDef(1).cellDef(1).initSpectralFlux = ones(1,p.nFineGroups);
% UO2 cell
g.regionDef(1).cellDef(2).name = 'UO2';
g.regionDef(1).cellDef(2).initZAIDs = [92235 92238 8016];
g.regionDef(1).cellDef(2).initNumDensities = [N_UO2(1), N_UO2(2), N_UO2(3)];
g.regionDef(1).cellDef(2).initSpectralFlux = ones(1,p.nFineGroups);

%% RUN VBUDSII
%cd(fullfile('..','..',''));
%addpath(fullfile('..','..',''));
addpath(p.vbudsiiDir);

[Results, p , g] = main(p,g);

r = Results.Region(1);

load xsprep1;

%% CHECK RESULTS
doPlot = 0;

if sum(r.Cell(1).spectralFlux ~= phi(:,1)) || ...
    sum(r.Cell(2).spectralFlux ~= phi(:,2))
    disp('xsprep1 flux test FAILED');
    figure;
    subplot(1,2,1);
    loglog(p.fineGroupDef(1:end-1)'*[1 1],[r.Cell(1).spectralFlux ...
        phi(:,1)])
    title('H2O')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('vbudsii','vbudsi','Location','Best')
    subplot(1,2,2);
    loglog(p.fineGroupDef(1:end-1)'*[1 1],[r.Cell(2).spectralFlux ...
        phi(:,2)])
    title('UO2')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('vbudsii','vbudsi','Location','Best')
    [r.Cell(2).spectralFlux./phi(:,2) r.Cell(2).spectralFlux-phi(:,2)]
   % print('xsprep1_111024.eps','-depsc','-r300');
else
    disp('xsprep1 flux test passed');
%    [Region(1).Cell(1).spectralFlux == phi(:,1), ...
%    Region(1).Cell(2).spectralFlux == phi(:,2)]
end

%{
phi(:,1) = r.Cell(1).spectralFlux;
phi(:,2) = r.Cell(2).spectralFlux;
save xsprep1.mat phi;
%}
