% This tests Multicell only against Geoff's original code.
% must specify exactly what this tests: tests the modified Multicell(geoff's
% vbudsii) against Geoff's original VBUDSII. Not a test for physical accuracy.

% do some cd-funk to get the absolute path of the vbudsii installation, then cd
% back into the directory that this file is in.

% MULTICELL4 (110918)

doPrint = 0;

mycd = cd;
cd(fullfile('..','..',''));
vbudsiiDir = cd;
cd(mycd); % alternatively, mfilename

resolveXS = [1 1];
S0iterthresh = [.01 .00001];
for itr = 1:length(resolveXS)
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
           'resolveXS',resolveXS(itr),...
           'S0iterthresh',S0iterthresh(itr)); % this largely affects the validation effort.

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

% define
% one region with two cells: UO2 and H2O
g.nRegions = 1;
g.regionDef(1).name = 'campaign1';
g.regionDef(1).nCells = 2;
g.regionDef(1).relVolumes = [4-pi/4, pi/4];
% H2O cell
g.regionDef(1).cellDef(1).name = 'H2O';
g.regionDef(1).cellDef(1).initZAIDs = [222];
g.regionDef(1).cellDef(1).initNumDensities = [.145];
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

[Results, p , g] = main(p,g);

r = Results.Region(1);

%% CHECK RESULTS
doPlot = 0;
load sigma4.mat
for cellidx = 2
figure;
        loglog(p.fineGroupDef(1:end-1)'*[1 1 1 1 1], ...
            [r.Cell(cellidx).fine(p.myMT(7)).value ...
             sum(r.Cell(cellidx).fine(p.myMT(2)).value)' ...
             r.Cell(cellidx).fine(p.myMT(4)).value ...
             r.Cell(cellidx).fine(p.myMT(18)).value ...
             r.Cell(cellidx).fine(p.myMT(102)).value],'LineWidth',2)
hold on
        loglog(p.fineGroupDef(1:end-1)'*[1 1 1 1 1], ...
            [Sigma.T(cellidx,:)' Sigma.E(cellidx,:)' Sigma.I(cellidx,:)' ...
             Sigma.F(cellidx,:)' Sigma.G(cellidx,:)']);
        legend('andover t','andover e','andover i','andover f','andover g', ...
               'austin t', 'austin e', 'austin i', 'austin f', 'austin g',...
               'Location','Best');
        xlabel('E (eV)')
        ylabel('\Sigma (cm^{-1}')
        title(['Macroscopic cross sections for cell ' g.regionDef(1).cellDef(cellidx).name]);
        if doPrint == 1
        print(['../../../../doc/figures/s0iteration1_' num2str(itr) ...
        'macxs_110928.eps'],'-depsc','-r300');
        end
end

cellidx = 2;
MT = 18;
lowS0 = zeros(p.nFineGroups,1);
highS0 = zeros(p.nFineGroups,1);
for zaididx = 1:length(r.Cell(cellidx).ZAIDs)
    lowS0 = lowS0 + r.Cell(cellidx).fine(p.myMT(MT)).z(zaididx).t(:,1)* ...
        r.Cell(cellidx).numDensities(zaididx);
    highS0 = highS0 + r.Cell(cellidx).fine(p.myMT(MT)).z(zaididx).t(:,end)* ...
        r.Cell(cellidx).numDensities(zaididx);
end

figure;
loglog(p.fineGroupDef(1:end-1)'*[1 1],...
    [r.Cell(cellidx).fine(p.myMT(MT)).value lowS0],'LineWidth',3);
hold on;
loglog(p.fineGroupDef(1:end-1)',...
    [highS0]);
legend('interpolated','lowS0','highS0');
title('S0 bounds')
disp('STDDD')
std(r.Cell(cellidx).fine(p.myMT(102)).value - lowS0)


load pi4.mat
if PI ~= r.PI
    PI(:,:,40)
    r.PI(:,:,40)
    disp('multicell4 PI test FAILED');
else
    disp('multicell4 PI test passed');
end

for cellidx = 1:g.regionDef(1).nCells
    for zaididx = 1:length(r.Cell(cellidx).ZAIDs)
    figure;
    semilogx(p.fineGroupDef(1:end-1), r.Cell(cellidx).S0s(:,zaididx));
    hold on;
    semilogx(p.fineGroupDef([1 end]), [-1 -1]);
    semilogx(p.fineGroupDef([1 end]), [10 10]);
    axis([0 inf -2 11])
    title(['s0 plot for cell ' num2str(cellidx) '; zaididx ' num2str(zaididx)])
    end
end

load multicell4data.mat % uses same data as multicell3
phi = ans;
if sum(r.Cell(1).spectralFlux ~= phi(:,1)) || ...
    sum(r.Cell(2).spectralFlux ~= phi(:,2))
    disp('multicell4 flux test FAILED');
    figure;
    subplot(1,2,1);
    semilogx(p.fineGroupDef(1:end-1)'*[1 1],[r.Cell(1).spectralFlux phi(:,1)])
    title('H2O')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('andover','austin')
    subplot(1,2,2);
    semilogx(p.fineGroupDef(1:end-1)'*[1 1],[r.Cell(2).spectralFlux phi(:,2)])
    title('UO2')
    ylabel('\phi (# n cm^{-2} s^{-1})')
    xlabel('E (eV)')
    legend('andover','austin')
    r.Cell(2).spectralFlux./phi(:,2);
    if doPrint == 1
    print(['../../../../doc/figures/s0iteration1_' num2str(itr) ...
        'flux_110928.eps'],'-depsc','-r300');
    %print('multicell4_fluxcompare_110921.eps','-depsc','-r300');
    end
else
    disp('multicell4 flux test passed');
%    [Region(1).Cell(1).spectralFlux == phi(:,1), ...
%    Region(1).Cell(2).spectralFlux == phi(:,2)]
end

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

























