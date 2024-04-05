% s0investigation1

% 1 - try only one s0 value for everything and compare flux results to
% actually doing s0 interpolation. To do this I should really make sure
% that my makeLibraryTempFlag == 1 and makeLibraryTempFlag == 2 give the
% same cross sections for a choice of S0

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
           'resolveXS',0,...
           'immutableMyMTs',0,...
           'makeRealNuFission',0,...
           'constS0',-1,...
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

[Results] = main(p,g);

r1 = Results.Region(1);

p.constS0 = 10;
%p.makeLibraryTempFlag = 1;

[Results, p, g, L] = main(p,g);

r2 = Results.Region(1);

if r1.Cell(1).spectralFlux ~= r2.Cell(1).spectralFlux | ...
        r1.Cell(2).spectralFlux ~= r2.Cell(2).spectralFlux
    disp('s0investigation flux test FAILED');
    figure;
    subplot(1,2,1);
    loglog(p.fineGroupDef(1:end-1)'*[1 1],...
        [r1.Cell(1).spectralFlux r2.Cell(1).spectralFlux])
    title('H2O')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('s0 = -1','s0 = 10','Location','Best')
    subplot(1,2,2);
    loglog(p.fineGroupDef(1:end-1)'*[1 1],...
        [r1.Cell(2).spectralFlux r2.Cell(2).spectralFlux])
    title('UO2')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('s0 = -1','s0 = 10','Location','Best')
    er0 = [r1.Cell(1).spectralFlux./r2.Cell(1).spectralFlux ...
        r1.Cell(2).spectralFlux./r2.Cell(2).spectralFlux]
    print('s0investigation_fluxcompare_111026.eps','-depsc','-r300');
else
    disp('s0investigation flux test passed');
%    [Region(1).Cell(1).spectralFlux == phi(:,1), ...
%    Region(1).Cell(2).spectralFlux == phi(:,2)]
end

B1=   N_UO2(1)*L.z(L.ZAID(92235)).m(L.MT(102)).xs(:,L.T(600),L.S0(-1));
B2=   N_UO2(1)*L.z(L.ZAID(92235)).m(L.MT(102)).xs(:,L.T(600),L.S0(10));
A1=   N_UO2(2)*L.z(L.ZAID(92238)).m(L.MT(102)).xs(:,L.T(600),L.S0(-1));
A2=   N_UO2(2)*L.z(L.ZAID(92238)).m(L.MT(102)).xs(:,L.T(600),L.S0(10));
C1=   N_UO2(3)*L.z(L.ZAID(8016)).m(L.MT(102)).xs(:,L.T(600),L.S0(-1));
C2=   N_UO2(3)*L.z(L.ZAID(8016)).m(L.MT(102)).xs(:,L.T(600),L.S0(10));

plotxsn(0,1,'Cell 2 MT 102',...
        {'indep -1','vbudsii -1','indep 10','vbudsii 10'},'',...
        p.fineGroupDef,...
        A1+B1+C1,...
        r1.Cell(2).fine(L.MT(102)).value,...
        A2+B2+C2,...
        r2.Cell(2).fine(L.MT(102)).value);

print('s0investigation_macroxscompare_111026.eps','-depsc','-r300');

% reaction rates
for cellidx = 1:g.regionDef(1).nCells
    for mt = 1:L.MT
        fprintf('CELL %s MT  %d  RR %f',...
            g.regionDef(1).cellDef(cellidx).name,mt,...
            r1.Cell(cellidx).fine(L.MT(mt)).RR);
    end
end

% 2 - try limiting the number of s0 interpolating values to just -1 and 5
% and just -1 and 10, and see if we have enough data points for an OKAY
% interpolation.

    % is this a moot point if the flux does not change from changing s0?
