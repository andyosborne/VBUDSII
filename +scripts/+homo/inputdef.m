function [p, g] = inputdef(reactortype, homoset, varargin)

import scripts.homo.*;

    % ZAIDS
[thisdir thisfname thisext] = fileparts(mfilename('fullpath'));
GEOFFMATCARD = importdata(fullfile(thisdir, 'imfmats.txt'), ' ', 3);
GEOFFMATCARD = GEOFFMATCARD.data;

zaidswanted = [11023 94239 94240 94241 94242 40090 40094 39089]';
%zaidswanted = [8016 94239]; % 94240 94241 94242 40090 40094 39089]';
idxs = [];
for i = 1:length(zaidswanted)
    idxs = [idxs find(GEOFFMATCARD(:,1) == zaidswanted(i))];
end
GEOFFMATCARD = GEOFFMATCARD(idxs',:);

ZAIDsToFarm = GEOFFMATCARD(:,1)';
DENSITIES = GEOFFMATCARD(:,2);


%%

fineGroupDef = 10.^(-2.8:0.1:7); 
%Thermal 98g structure

%fineGroupDef = [501.187000000000,630.957000000000,794.328000000000,1000,1258.93000000000,1584.89000000000,1995.26000000000,2511.89000000000,3162.28000000000,3981.07000000000,5011.87000000000,6309.57000000000,7943.28000000000,10000,12589.3000000000,15848.9000000000,19952.6000000000,25118.9000000000,31622.8000000000,39810.7000000000,50118.7000000000,63095.7000000000,79432.8000000000,100000,125893,158489,199526,251189,316228,398107,501187,630957,794328,1000000,1258930,1584890,1995260,2511890,3162280,3981070,5011870,6309570,7943280 8e6 9e6 10e6 11e6 12e6 13e6 14e6 15e6 16e6 17e6 18e6 19e6 20e6];
%Fast 56g structure

%fineGroupDef=[501.187000000000,630.957000000000,794.328000000000,1000,1258.93000000000,1584.89000000000,1995.26000000000,2511.89000000000,3162.28000000000,3981.07000000000,5011.87000000000,6309.57000000000,7943.28000000000,10000,12589.3000000000,15848.9000000000,19952.6000000000,25118.9000000000,31622.8000000000,39810.7000000000,50118.7000000000,63095.7000000000,79432.8000000000,100000,125893,158489,199526,251189,316228,398107,501187,630957,794328,1000000,1258930,1584890,1995260,2511890,3162280,3981070,5011870,6309570,7943280];
%Fast 42g structure


nFineGroups = length(fineGroupDef) - 1;




%%
% These group bounds are provided by the Multi benchmark.
fewGroupDef = fineGroupDef([1 ...
    find(fineGroupDef > 4, 1, 'first'), ...
    find(fineGroupDef > 9.1e3, 1, 'first') ...
    end]);
%Thermal fine group structure



% fewGroupDef = fineGroupDef([1 ...
%     find(fineGroupDef > 10000, 1, 'first'), ...
%     find(fineGroupDef > 251189, 1, 'first') ...
%     end]);
%Fast Fine Group structure

%%


% VBUDSII parameters.
p = vbudsii.initparam();
p.nFineGroups = nFineGroups;
p.nFewGroups = length(fewGroupDef) - 1;
p.fineGroupDef = fineGroupDef;
p.fewGroupDef = fewGroupDef;
p.XSLibraryMAT = 'XSLibraryHomo.mat';
p.makeLibrary = 0;
%HELLO FUTURE ME (10.5.18), LETS FIND OUT WHAT MAKELIBRARY DOES
%makeLibrary runs that download function (But only sorta)
%p.makeLibraryTempFlag = makelib([ZAIDsToFarm 24052 26056 28058], fineGroupDef);
%HELLO FUTURE ME. THIS WAS COMMENTED OUT IN THE ORIGINAL FOR SOME REASON
p.immutableMyMTs = 2; % Decides when derived xs's are made.
p.useWatt = 1;
p.totalfortransport = 1;
p.stammlertransport = 1;
p.resolveXS = 1;

if strcmp(homoset, 'off')
elseif strcmp(homoset, 'orig')
    p.homogenize = 1;
    p.homogenizeThreshFac = .8; %threshfac;
    p.homogenizeSlopeFac = .5; %slopefac;
elseif strcmp(homoset, 'indiv')
    p.homogenize = 2;
    p.homogenizeThreshFac = .8; %threshfac;
    p.homogenizeSlopeFac = .5; %slopefac;
else
    error('Cooperate with me here! Not much time left.');
    %lol
end

% VBUDSII geometry.
g = vbudsii.initgeom();

uc = struct();


if strcmp(reactortype, 'fastgeoff')
    % Determining pin diameter.
    outrad = 0.538;
    effPitch = sqrt(pi * outrad^2);
    uc.pinPitch = effPitch * [1 1]; % % cm

    uc.pinDiam = 2*0.289*[1 1]; %sqrt(uc.pinPitch^2 * fuelVolFraction * 4 / pi);
else
    % Determining pin diameter.
    uc.pinPitch = [1.2 1.2]; % % cm

    uc.pinDiam = [0.8 0.8]; %sqrt(uc.pinPitch^2 * fuelVolFraction * 4 / pi);
    rodDiam = 0.9;
end

uc.f = [0.75 0.25]; % some weighting of the fuel regions
uc.g = [0.75 0.25]; % some weighting of the moderator regions
uc.sauerConst.mod = 2.35;
uc.sauerConst.fuel = 5.00;

% define
% one region with two cells: UO2 and H2O
g.nRegions = 1;
g.regionDef(1).name = 'campaign1';
g.regionDef(1).nCells = 4;
g.regionDef(1).relVolumes = [...
    pi/4*uc.pinDiam(1)^2,...
    uc.pinPitch(1)^2-pi/4*uc.pinDiam(1)^2 ...
    pi/4*uc.pinDiam(2)^2,...
    uc.pinPitch(2)^2-pi/4*uc.pinDiam(2)^2 ...
    ];

% possible materials.
% natural uo2
density_UO2 = 11;
[ao, wo, N_UO2] = matl([92235 0.0072;
                        92238 1-.0072;
                        8016  2], 1, density_UO2);
nuo2cell.name = '';
nuo2cell.isFissionable = 0;
nuo2cell.initZAIDs = [92235 92238 8016];
nuo2cell.initNumDensities = N_UO2';
nuo2cell.initDensity = density_UO2;
nuo2cell.initTemp = 294;

% enriched uo2
UO2_composition = [92235 0.00070803;
                   92238 0.022604;
                   8016  0.046624];
[ao, wo, N_UO2, density_UO2] = matl(UO2_composition, 3);
euo2cell.name = '';
euo2cell.isFissionable = 0;
euo2cell.initZAIDs = [92235 92238 8016];
euo2cell.initNumDensities = N_UO2';
euo2cell.initDensity = density_UO2;
euo2cell.initTemp = 294;

% highly enriched uo2
UO2_composition = [92235 8 * 0.00070803;
                   92238 0.022604;
                   8016  0.046624];
[ao, wo, N_UO2, density_UO2] = matl(UO2_composition, 3);
heucell.name = '';
heucell.isFissionable = 0;
heucell.initZAIDs = [92235 92238 8016];
heucell.initNumDensities = N_UO2';
heucell.initDensity = density_UO2;
heucell.initTemp = 294;

% mildly enriched uo2
UO2_composition = [92235 4 * 0.00070803;
                   92238 0.022604;
                   8016  0.046624];
[ao, wo, N_UO2, density_UO2] = matl(UO2_composition, 3);
meucell.name = '';
meucell.isFissionable = 0;
meucell.initZAIDs = [92235 92238 8016];
meucell.initNumDensities = N_UO2';
meucell.initDensity = density_UO2;
meucell.initTemp = 294;

% pure uo2
UO2_composition = [92238 0.022604+0.00070803;
                   8016  0.046624];
[ao, wo, N_UO2, density_UO2] = matl(UO2_composition, 3);
puo2cell.name = '';
puo2cell.isFissionable = 0;
puo2cell.initZAIDs = [92238 8016];
puo2cell.initNumDensities = N_UO2';
puo2cell.initDensity = density_UO2;
puo2cell.initTemp = 294;

% IMF
density_IMF = 5.6880202; %g/cm^3
[ao, wo, N_IMF, densities_IMF] = matl(GEOFFMATCARD, 2, density_IMF);
imfcell.name = '';
imfcell.isFissionable = 0;
imfcell.initZAIDs = ZAIDsToFarm;
imfcell.initNumDensities = N_IMF';
imfcell.initDensity = density_IMF;
imfcell.initTemp = 294;

% pure PuO2
density_PuO2 = 11.5;
[ao, wo, N_PuO2, densities_PuO2] = matl([94239 1;
                                         8016 2], 1, density_PuO2);
plutcell.name = '';
plutcell.isFissionable = 0;
plutcell.initZAIDs = [94239 8016];
plutcell.initNumDensities = [N_PuO2'];
plutcell.initDensity = density_PuO2;
plutcell.initTemp = 294;

% geoff puo2
density_PuO2 = 11.5;
[ao, wo, N_PuO2, densities_PuO2] = matl([94239 .5;
                                         94240 .5;
                                         8016 2], 1, density_PuO2);
                                         %94241 .2;
plutgcell.name = '';
plutgcell.isFissionable = 0;
%plutgcell.initZAIDs = [94239 94240 94241 8016];
plutgcell.initZAIDs = [94239 94240 8016];
plutgcell.initNumDensities = [N_PuO2'];
plutgcell.initDensity = density_PuO2;
plutgcell.initTemp = 294;

% water
H2O_composition = [1001 0.066988;
                   8016 0.033414];
[ao, wo, N_H2O, density_H2O] = matl(H2O_composition, 3);
basicwatercell.name = '';
basicwatercell.isFissionable = 0;
basicwatercell.initZAIDs = [11 8016];
basicwatercell.initNumDensities = [N_H2O];
basicwatercell.initDensity = density_H2O;
basicwatercell.initTemp = 294;

% sodium
density_Na = 0.882; % g/cm^3
[ao, wo, N_Na] = matl([11023 1], 1, density_Na);
sodiumcell.name = '';
sodiumcell.isFissionable = 0;
sodiumcell.initZAIDs = [11023];
sodiumcell.initNumDensities = [N_Na'];
sodiumcell.initDensity = density_Na;
sodiumcell.initTemp = 294;

% sodiumsmear
[ao, wo, N_Annulus, density_Ann] = ...
    matl([11023 1.038e-2;
          24052 4.973e-3;
          26056 1.728e-2;
          28058 3.627e-3
        ], 3);
sodiumsmearcell.name = '';
sodiumsmearcell.isFissionable = 0;
sodiumsmearcell.initZAIDs = [11023 24052 26056 28058];
sodiumsmearcell.initNumDensities = [N_Annulus'];
sodiumsmearcell.initDensity = density_Ann;
sodiumsmearcell.initTemp = 294;


if strcmp(reactortype, 'equi')

    g.regionDef(1).cellDef(1) = euo2cell;
    g.regionDef(1).cellDef(1).name = 'UOX1';
    g.regionDef(1).cellDef(3) = euo2cell;
    g.regionDef(1).cellDef(3).name = 'UOX2';

    g.regionDef(1).cellDef(2) = basicwatercell;
    g.regionDef(1).cellDef(2).name ='UOX1water';
    g.regionDef(1).cellDef(4) = basicwatercell;
    g.regionDef(1).cellDef(4).name ='UOX2water';

    uc.PinCellNames = {'UOX1', 'UOX2'};
    uc.AntipinCellNames = {'UOX1water', 'UOX2water'};

elseif strcmp(reactortype, 'pure')
    
    g.regionDef(1).cellDef(1) = euo2cell;
    g.regionDef(1).cellDef(1).name = 'UOX';
    g.regionDef(1).cellDef(3) = puo2cell;
    g.regionDef(1).cellDef(3).name = 'PureUOX';

    g.regionDef(1).cellDef(2) = basicwatercell;
    g.regionDef(1).cellDef(2).name ='UOXwater';
    g.regionDef(1).cellDef(4) = basicwatercell;
    g.regionDef(1).cellDef(4).name ='PureUOXwater';

    uc.PinCellNames = {'UOX', 'PureUOX'};
    uc.AntipinCellNames = {'UOXwater', 'PureUOXwater'};

elseif strcmp(reactortype, 'imf')

    g.regionDef(1).cellDef(1) = euo2cell;
    g.regionDef(1).cellDef(1).name = 'UOX';
    g.regionDef(1).cellDef(3) = imfcell;
    g.regionDef(1).cellDef(3).name = 'IMF';

    g.regionDef(1).cellDef(2) = basicwatercell;
    g.regionDef(1).cellDef(2).name ='UOX water';
    g.regionDef(1).cellDef(4) = basicwatercell;
    g.regionDef(1).cellDef(4).name ='IMF water';

    uc.PinCellNames = {'UOX', 'IMF'};
    uc.AntipinCellNames = {'UOX water', 'IMF water'};

elseif strcmp(reactortype, 'fast')

    g.regionDef(1).cellDef(1) = heucell;
    g.regionDef(1).cellDef(1).name = 'HEU';
    g.regionDef(1).cellDef(3) = meucell;
    g.regionDef(1).cellDef(3).name = 'MEU';

    g.regionDef(1).cellDef(2) = sodiumcell;
    g.regionDef(1).cellDef(2).name = 'HEUsodium';
    g.regionDef(1).cellDef(4) = sodiumcell;
    g.regionDef(1).cellDef(4).name = 'MEUsodium';

    uc.PinCellNames = {'HEU', 'MEU'};
    uc.AntipinCellNames = {'HEUsodium', 'MEUsodium'};

elseif strcmp(reactortype, 'fastplut')

    g.regionDef(1).cellDef(1) = heucell;
    g.regionDef(1).cellDef(1).name = 'HEU';
    g.regionDef(1).cellDef(3) = plutcell;
    g.regionDef(1).cellDef(3).name = 'PuO$_2$';

    g.regionDef(1).cellDef(2) = sodiumcell;
    g.regionDef(1).cellDef(2).name = 'HEU sodium';
    g.regionDef(1).cellDef(4) = sodiumcell;
    g.regionDef(1).cellDef(4).name = 'PuO$_2$ sodium';

    uc.PinCellNames = {'HEU', 'PuO$_2$'};
    uc.AntipinCellNames = {'HEU sodium', 'PuO$_2$ sodium'};

elseif strcmp(reactortype, 'fastgeoff')

    g.regionDef(1).cellDef(1) = nuo2cell;
    g.regionDef(1).cellDef(1).name = 'UOX fuel';
    g.regionDef(1).cellDef(3) = plutgcell;
    g.regionDef(1).cellDef(3).name = 'PuO$_2$ fuel';

    g.regionDef(1).cellDef(2) = sodiumsmearcell;
    g.regionDef(1).cellDef(2).name = 'UOX coolant';
    g.regionDef(1).cellDef(4) = sodiumsmearcell;
    g.regionDef(1).cellDef(4).name = 'PuO$_2$ coolant';

    uc.PinCellNames = {'UOX fuel', 'PuO$_2$ fuel'};
    uc.AntipinCellNames = {'UOX coolant', 'PuO$_2$ coolant'};


else
    error('Who do you think you are?');
    %LOL
    %I can tell most of this was written at 3am
end

g.regionDef(1).uc = uc;











return;















% UO2 cell
g.regionDef(1).cellDef(1).name = 'UOX fuel';
g.regionDef(1).cellDef(1).initZAIDs = [92235 92238 8016];
g.regionDef(1).cellDef(1).initNumDensities = N_UO2';
g.regionDef(1).cellDef(1).initDensity = density_UO2;
g.regionDef(1).cellDef(1).initTemp = 294;

% IMF cell
g.regionDef(1).cellDef(3).name = 'IMF fuel';
if ~useUO2
    g.regionDef(1).cellDef(3).initZAIDs = ZAIDsToFarm;
    g.regionDef(1).cellDef(3).initNumDensities = N_IMF';
    g.regionDef(1).cellDef(3).initDensity = density_IMF;
    g.regionDef(1).cellDef(3).initTemp = 294;
else
    UO2_composition = [92235 ratio*0.00070803;
    92238 0.022604;
    8016  0.046624];
    [ao, wo, N_UO2, density_UO2] = matl(UO2_composition, 3);
    g.regionDef(1).cellDef(3).initZAIDs = [92235 92238 8016];
    g.regionDef(1).cellDef(3).initNumDensities = N_UO2';
    g.regionDef(1).cellDef(3).initDensity = density_UO2;
    g.regionDef(1).cellDef(3).initTemp = 294;
end
Z_composition = [40090 0.043241];

g.regionDef(1).cellDef(2).name = 'UOX coolant';
g.regionDef(1).cellDef(4).name = 'IMF coolant';

    % Just ignore the zirconium.
    [ao, wo, N_H2O, density_H2O] = matl(H2O_composition, 3);
    g.regionDef(1).cellDef(2).initZAIDs = [11 8016];
    g.regionDef(1).cellDef(2).initNumDensities = [N_H2O];
    g.regionDef(1).cellDef(2).initDensity = density_H2O;
    g.regionDef(1).cellDef(2).initTemp = 294;

    g.regionDef(1).cellDef(4).initZAIDs = [11 8016];
    g.regionDef(1).cellDef(4).initNumDensities = [N_H2O];
    g.regionDef(1).cellDef(4).initDensity = density_H2O;
    g.regionDef(1).cellDef(4).initTemp = 294;















% FarmXS input.

%s0values = [10 3.5 1.5 1 .5 0 -.5 -1];
%for idx = 1:length(s0values)
%    s0values(idx) = str2double(sprintf('%.1e', 10^(s0values(idx))));
%end
s0values = 10.^[10 5 3 2 1 0 -1];
pL = struct(...
    'metamat', 'FarmXSMetaHomo4', ...
    'address', fullfile(thisdir, 'XSFarm'), ...
    'verbose', 1, ...
    'explore', 1, ...
'endfversion', 7.1, ...
       'mine', 1, ...
      'plant', 1, ...
  'fertilize', 1, ...
    'harvest', 1, ...
          'e', struct('searchWeb', 0), ...
          'm', struct('downloadTape20s', 1), ...
          'p', struct('MFMTs', [3 1;
                                3 18;
                                3 102;
                                3 251;
                                3 452;
                                ], ...
                        'MFs', [6], ...
                         'Ts', [294], ...
                        'S0s', s0values, ...
                     'IGNstr', '1', ... % read in group structure
                   'groupDef', fineGroupDef, ...
                     'IWTstr', '3'), ... % maybe 7 for fast
           'f', struct('overwriteOutput', 1, ...
                            'logResults', 1), ...
           'h', struct('makeLibrary', 1, ...
                        'makeArrays', 0), ...
                 'ZAIDs',[ZAIDsToFarm 11 1001 8016 92235 92238]);
       %'ZAIDs', [11 1001 8016 40090 92235 92238]);




