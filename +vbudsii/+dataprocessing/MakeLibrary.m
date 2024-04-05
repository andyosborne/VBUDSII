function L = MakeLibrary(p)
%
%VBUDSII/DATAPROCESSING/MAKELIBRARY Initializes and populates the "static" data
% structure L, which contains data pulled from ENDF/B via NJOY.
%
% INPUT
%   p           Parameter structure.
%
% OUTPUT
%   L           Cross section/data library containing data from ENDF/B and
%               processed by NJOY. The structure is initially created in the
%               MakeLibrary.m function. The script that calls VBUDSII may also
%               contain code that analyzes Results, and accordingly it may be
%               useful to have the Library available for that analysis.
%
% DEPENDENCIES
%   VBUDSII/DATA/FULLLIBRARY.mat
%   VBUDSII/DATA/water.mat
%   VBUDSII/DATA/reducedlib.mat

% NOTES
%   The current state of this function is such that it is adapting to a cross
%   section grabbing process that is not complete. It adapts to a method that
%   allows verification against the previous version of VBUDSI, allows
%   connecting directly to a library created via NJOY scripts, and also
%   connects to pre-made MAT-files that contain cross sections under a
%   different naming convention. The ability to verify VBUDSII against other
%   codes depends on making sure that the Library is populated as desired. It
%   is possible that this entire function can be done away with, once the
%   correct process for obtaining cross sections is developed.
%
% MAJOR REVISIONS
%   date        handle      description
%   20111030    cld2469     writing comments
%
% TASKLIST
%   1- Finish generating cross sections properly through NJOY scripts.
%   2- Make permanent the immutable MTs 6 - 9.
%   3- Water's transport cross section is wrong so long as the 1001 used for it
%   does not come from H1w.

import vbudsii.*

util.PrintEntering(p, 'MakeLibrary');
    
% Get the directory this file is in.
S = dbstack('-completenames');
thisdir = fileparts(S(1).file);

if isstruct(p.makeLibraryTempFlag)

    if isfield(p.makeLibraryTempFlag, 'endfversion') && ...
            p.makeLibraryTempFlag.endfversion == 7.1
        disp(p.makeLibraryTempFlag);
        [Meta, L] = FarmXSVIIp1(p.makeLibraryTempFlag);
    else
        disp(p.makeLibraryTempFlag)
        [Meta, L] = FarmXS(p.makeLibraryTempFlag);
    end

    % create water. eventually want to move away from this dependency.
    % requires that I use Geoff's temperatures, unless I also provide a Tidxs.
    load(fullfile(thisdir, '..', '+data', ...
        'XSLibrary_water_byGeoff_120101.mat'));

    S0idxs = [1:5 7 12];
    Tidxs = [1 4 6 8]; % 296 350 400 450 500 600 800 1000
    MCNPXcorrectionFactor = 1; %%% 1.888
    if isfield(p,'MCNPXwaterCorrectionFactor')
        MCNPXcorrectionFactor = p.MCNPXwaterCorrectionFactor;
        disp('NOTE: Using MCNPXwaterCorrectionFactor.');
    end
    %for mtidx = [3 5] % these two indices pick out MT 102 and MT 251
    % Geoff has added the inelastic cross section to MT 2 in his parse_bound.
    % 1 2 102 222 251
    % 251 seems to be all zero.
    if isfield(p, 'MCNPXcorrectTheseMTIDXS')
        mtidxsToCorrect = p.MCNPXcorrectTheseMTIDXS;
        disp('NOTE: MCNPXcorrectTheseMTIDXS');
    else
        mtidxsToCorrect = [3 5];
    end
    for mtidx = [3 5] % these two indices pick out MT 102 and MT 251
        if sum(mtidx == mtidxsToCorrect)
            vectorscale = MCNPXcorrectionFactor;
        else
            vectorscale = 1;
        end
        if L.ZAIDs == 222
            L.z(L.ZAID(222)).m(L.MT(light{1}.MT(mtidx))).xs = ...
                 vectorscale * light{5}.m(mtidx).xs(:,Tidxs,S0idxs);
        end
    end
    if (isfield(p, 'MCNPXcorrectKernel') && p.MCNPXcorrectKernel) || ...
            ~isfield(p, 'MCNPXcorrectKernel')
        kernelscale = MCNPXcorrectionFactor;
    else
        kernelscale = 1;
    end

    if L.ZAIDs == 222
        % Scale the kernel.
        L.z(L.ZAID(222)).m(L.MT(2)).xs = kernelscale * ...
                light{5}.m(2).kern(:,:,Tidxs,S0idxs);

    end

    if isfield(p, 'waterPiece') && p.waterPiece
        load(fullfile(thisdir, '..', '+data', ...
            'lwtrscatteringkernel120624'));
        L.z(L.ZAID(1001)).m(L.MT(2)).xs = ...
            lwtrscattering(:,:, Tidxs, S0idxs);
    end

    if isfield(p, 'waterPieceEquiprob') && p.waterPieceEquiprob
        load(fullfile(thisdir, '..', '+data', ...
            'lwtrscatteringequiprob'));
        L.z(L.ZAID(1001)).m(L.MT(2)).xs = ...
            equiprob(:,:, Tidxs, S0idxs);
    end

    if isfield(p, 'waterWhole') && p.waterWhole
        load(fullfile(thisdir, '..', '+data', ...
            'waternuclide120624'));
        L.z(L.ZAID(222)).m(L.MT(2)).xs = ...
            waterscatteringkernel(:,:,Tidxs,S0idxs);
        L.z(L.ZAID(222)).m(L.MT(102)).xs = ...
            waterabsorption(:,Tidxs,S0idxs);
    end

    % Create an inelastic kernel from the sum of its parts.
    if isfield(p, 'noInelastic')
        disp('NOTE: USING noInelastic.');
    else
        for z = L.ZAIDs
            if z ~= 222
                inelasticSum = 0;
                for m = L.z(L.ZAID(z)).inelasticMTs
                    L.z(L.ZAID(z)).m(L.MT(6)).xs = L.z(L.ZAID(z)).m(L.MT(6)).xs + ...
                        L.z(L.ZAID(z)).m(L.MT(m)).xs;
                end
                L.z(L.ZAID(z)).m(L.MT(2)).xs = L.z(L.ZAID(z)).m(L.MT(2)).xs + ...
                    L.z(L.ZAID(z)).m(L.MT(6)).xs;
            end
        end
    end


    % L.MTs = [1 2 4 6:9 18 102 251 452]; depends on how this works in
    % FarmXs

end

%% Precompute my MTs 6 - 9. This should eventually be permanent.
% We should do this here to avoid making a "macroscopic" nu, or mubar.
if p.immutableMyMTs == 1

    if isfield(p, 'zeromu') && p.zeromu
        L.z(L.ZAID(1001)).m(L.MT(251)).xs(1:26,:,:) = ...
            zeros(26,length(L.Ts), length(L.S0s));
    end

    % For each ZAID.
    for z = L.ZAIDs

        % absorption
        absorptionxs = ...
            L.z(L.ZAID(z)).m(L.MT(18)).xs + L.z(L.ZAID(z)).m(L.MT(102)).xs;

        % MT = 9, nu-fission
        if p.makeRealNuFission % && L.z(L.ZAID(z)).isFissionable
            % Actually use nu.

            L.z(L.ZAID(z)).m(L.MT(9)).xs = ...
                  L.z(L.ZAID(z)).m(L.MT(452)).xs .* ...
                  L.z(L.ZAID(z)).m(L.MT(18)).xs;

        else
            % Don't use nu -- to find agreement with Geoff's VBUDSII.
            L.z(L.ZAID(z)).m(L.MT(9)).xs = L.z(L.ZAID(z)).m(L.MT(18)).xs;

        end

        for t = L.Ts
            for s0 = L.S0s
                % MT = 7, total.
                L.z(L.ZAID(z)).m(L.MT(7)).xs(:,L.T(t),L.S0(s0)) = ...
                    sum(L.z(L.ZAID(z)).m(L.MT(2)).xs(:,:,L.T(t),L.S0(s0)))' ...
                    + absorptionxs(:,L.T(t),L.S0(s0));
                 %  + L.z(L.ZAID(z)).m(L.MT(4)).xs(:,L.T(t),L.S0(s0))...
                 %  + (L.z(L.ZAID(z)).m(L.MT(4)).xs(:,L.T(t),L.S0(s0)) - ...
                 %  sum(L.z(L.ZAID(z)).m(L.MT(6)).xs(:,:,L.T(t),L.S0(s0)))');
%                    + L.z(L.ZAID(z)).m(L.MT(6)).xs(:,L.T(t),L.S0(s0));
%

                if z == 222 && ((isfield(p,'waterTransportFromParts') && ...
                        p.waterTransportFromParts == 1) || ...
                        ~isfield(p,'waterTransportFromParts'))
                else
                    % Do not make transport like this for water.
                    % MT = 8, transport.
                    L.z(L.ZAID(z)).m(L.MT(8)).xs(:,L.T(t),L.S0(s0)) = ...
                        L.z(L.ZAID(z)).m(L.MT(7)).xs(:,L.T(t),L.S0(s0)) - ...
                        L.z(L.ZAID(z)).m(L.MT(251)).xs(:,L.T(t),L.S0(s0)) .* ...
                        sum(L.z(L.ZAID(z)).m(L.MT(2)).xs(:,:,L.T(t),L.S0(s0)))';
                end

            end
        end
    end

    % Create water's transport cross section "from parts".
    if sum(L.ZAIDs == 222) && ((isfield(p,'waterTransportFromParts') && ...
            p.waterTransportFromParts == 1) || ...
            ~isfield(p,'waterTransportFromParts'))
        L.z(L.ZAID(222)).m(L.MT(8)).xs = 2*L.z(L.ZAID(1001)).m(L.MT(8)).xs + ...
                                           L.z(L.ZAID(8016)).m(L.MT(8)).xs;
    end
end

%% S0 DEPENDENCE
%% Print a statement to the command window for each ZAID-MT that shows
%% dependency on S0, given the way that cross sections are currently generated.
%fprintf('\nCHECKING S0 SERIOUSNESS: THESE XS''S DEPEND ON S0\n');
%for z = L.ZAIDs
%    for m = L.mainMTs
%        for t = L.Ts
%            % Do not look at scattering.
%            if (m ~= 2) && ...
%                    ( sum(abs(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,L.T(t),1) ...
%                     -L.z(L.ZAID(z)).m(L.MT(m)).xs(:,L.T(t),end))) > 0 )
%                 % S0 = -1 to S0 = 10 yields different cross section values.
%
%                 % Print a statement to the command window.
%                 fprintf('ZAID %d   MT %d    T %d \n',z,m,t)
%            end
%        end
%    end
%end

util.PrintExiting(p, 'MakeLibrary');

end



function L = initL(p)
% Library.z(L.ZAID(222)).m(L.MT(18)).xs(energy,L.T(300),L.S0(-1))

    ZAIDs = [1001 8016 11023 92235 92238 222];
    MTs = [2 4 6 7 8 9 18 102 251 452];
    Ts = [300 600 900 1200 1500];
    S0s = [-1 0 1 2 3 5 10];

    L = struct('ZAIDs',ZAIDs,...
               'MTs',MTs,...
               'Ts',Ts,...
               'S0s',S0s,...
               'ZAID',[],...
               'MT',[],...
               'T',[],...
               'S0',[],...
        'groupDef',p.fineGroupDef,...
               'z',struct('isFissionable',[],...
                          'm',struct('hasResonances',[],...
                                     'xs',[])));

    L.ZAID(ZAIDs) = 1:length(ZAIDs);
    L.MT(MTs) = 1:length(MTs);
    L.T(Ts) = 1:length(Ts);
    L.S0 = containers.Map({-1 0 1 2 3 5 10},{1 2 3 4 5 6 7});

    L.ZAID = sparse(L.ZAID);
    L.MT = sparse(L.MT);
    L.T = sparse(L.T);
    nTs = length(Ts);
    nS0s = length(S0s);

    for i = ZAIDs
        L.z(L.ZAID(i)).isFissionable = 0;
        for j = MTs
            L.z(L.ZAID(i)).m(L.MT(j)).hasResonances = 1;
            if j == 2 || j == 6
                L.z(L.ZAID(i)).m(L.MT(j)).xs = ...
                    zeros(p.nFineGroups,p.nFineGroups,nTs,nS0s);
            %elseif j == 452
            %% THIS IS TEMPORARY, SINCE ANDY DOES NOT STORE 452
            %    L.z(L.ZAID(i)).m(L.MT(j)).xs = ones(p.nFineGroups,nTs,nS0s);
            else
                L.z(L.ZAID(i)).m(L.MT(j)).xs = zeros(p.nFineGroups,nTs,nS0s);
            end
        end
    end

end

function L = initLibrary(p)

    ZAIDs = [1001 8016 11023 92235 92238 222];
    Ts = [300 600 900 1200 1500];
    S0s = [-1 0 1 2 3 5 10];

mainMTs = [1 2 6 7 8 9 18 102 251 452];
MTs = mainMTs;

L = struct('groupDef',p.fineGroupDef,...
    'nGroups',length(p.fineGroupDef)-1,...
    'ZAIDs',ZAIDs,...
    'MTs',MTs,...
    'mainMTs',mainMTs,...
    'Ts',Ts,...
    'S0s',S0s,...
    'ZAID',[],...
    'MT',[],...
    'mainMT',[],...
    'T',[],...
    'S0',[],...
    'MTmap',[],...
    'z',struct('isFissionable',[],...
    'inelasticMTs',[],...
    'm',struct('hasResonances',[],...
    'xs',[])));
%           'MTnames', cell(1,length(MTs)),...

L.ZAID(ZAIDs) = 1:length(ZAIDs);
L.MT(MTs) = 1:length(MTs);
L.mainMT(MTs) = 1:length(MTs);
L.T(Ts) = 1:length(Ts);
L.S0 = containers.Map({-1 0 1 2 3 5 10},{1 2 3 4 5 6 7});

L.ZAID = sparse(L.ZAID);
L.MT = sparse(L.MT);
L.mainMT = sparse(L.mainMT);
L.T = sparse(L.T);
nTs = length(Ts);
nS0s = length(S0s);

for i = ZAIDs
    L.z(L.ZAID(i)).isFissionable = 0;
    for j = MTs
        L.z(L.ZAID(i)).m(L.MT(j)).hasResonances = 1;
        if j == 2  || j == 6 || j == 221 || (j >= 51 && j <= 91) % VERY DANGEROUS. EDIT.
            L.z(L.ZAID(i)).m(L.MT(j)).xs = ...
                zeros(L.nGroups,L.nGroups,nTs,nS0s);
        %elseif j == 452
        %    % THIS IS TEMPORARY, SINCE ANDY DOES NOT STORE 452
        %    L.z(L.ZAID(i)).m(L.MT(j)).xs = ones(L.nGroups,nTs,nS0s);
        else
            L.z(L.ZAID(i)).m(L.MT(j)).xs = zeros(L.nGroups,nTs,nS0s);
        end
    end
end

MTcellmap1 = {1;
    2;
    4;
    6;
    7;
    8;
    9;
    18;
    102;
    251;
    452};
MTcellmap2 = {'total';
    'elastickernel';
    'inelastic';
    'vbudsiiinelastickernel';
    'vbudsiitotal';
    'vbudsiitransport';
    'vbudsiinufission';
    'fission';
    'gamma';
    'mubar (average scattering cosine)';
    'nu (neutrons per fission)'};

L.MTmap = containers.Map(MTcellmap1, MTcellmap2);

end















































%{
lib = CrossSectionLibrary.soleInstance();

if ~isstruct(p.makeLibraryTempFlag) && p.makeLibraryTempFlag == 0
    % This option is the ideal usage, where the Library is loaded from NJOY
    % scripts. However, hopefully water would not need to be added like it is
    % here.

    load(p.XSLibraryMAT);

    % Just add water. First, account for this with the indexers.
    L.ZAIDs = [L.ZAIDs 222];
    L.ZAID(L.ZAIDs) = 1:length(L.ZAIDs);
    nTs = length(L.Ts);
    nSOs = length(L.S0s);

    % Initialize new data fields.
    L.z(L.ZAID(222)).m(L.MT(2)).xs = ...
        zeros(p.nFineGroups,p.nFineGroups,nTs,nSOs);

    L.z(L.ZAID(222)).m(L.MT(4)).xs = zeros(p.nFineGroups,nTs,nSOs);
    L.z(L.ZAID(222)).m(L.MT(18)).xs = zeros(p.nFineGroups,nTs,nSOs);
    L.z(L.ZAID(222)).m(L.MT(102)).xs = zeros(p.nFineGroups,nTs,nSOs);

    L.z(L.ZAID(222)).m(L.MT(251)).xs = zeros(p.nFineGroups,nTs,nSOs);

    L.z(L.ZAID(222)).isFissionable = 0;

    % Assign cross sections to data fields. There is no S0 resolution here.
    for temp = L.Ts
        for s0 = L.S0s
            L.z(L.ZAID(222)).m(L.MT(2)).xs(:,:,L.T(temp),L.S0(s0)) = ...
                lib.elScatKernel('H2O', temp);

            L.z(L.ZAID(222)).m(L.MT(4)).xs(:,L.T(temp),L.S0(s0)) = ...
                lib.xsSpectrum('H2O', 'inelastic', temp);
            L.z(L.ZAID(222)).m(L.MT(18)).xs(:,L.T(temp),L.S0(s0)) = ...
                lib.xsSpectrum('H2O', 'fission', temp);
            L.z(L.ZAID(222)).m(L.MT(102)).xs(:,L.T(temp),L.S0(s0)) = ...
                lib.xsSpectrum('H2O', 'nGamma', temp);

            L.z(L.ZAID(222)).m(L.MT(251)).xs(:,L.T(temp),L.S0(s0)) = ...
                lib.muScat('H2O', temp);
        end
    end

elseif ~isstruct(p.makeLibraryTempFlag) && p.makeLibraryTempFlag == 1
    % This is legacy, and is really only used to compare to Geoff's VBUDSII
    % results. The cross sections are populated from Andy's OO cross section
    % code,

    % Subfunction below.
    L = initLibrary(p);

    % Create cell arrays for indexing purposes, between the old XS format and
    % the new one.
    z = {1001, 'H1';
         8016, 'O16';
         92235, 'U235';
         92238, 'U238';
         222, 'H2O'};

    m = {4,'inelastic';
         18,'fission';
         102,'nGamma'};

    % For each ZAID.
    for i = 1:5

        if z{i,1} == 92235 || z{i,1} == 92238
            L.z(L.ZAID(z{i,1})).isFissionable = 1;
        end

        for t = L.Ts
            for b = L.S0s

                % For each MT.
                for j = 1:3
                  L.z(L.ZAID(z{i,1})).m(L.MT(m{j,1})).xs(:,L.T(t),L.S0(b)) ...
                    = lib.xsSpectrum(z{i,2},m{j,2},t);
                end

                % Unique data.
                L.z(L.ZAID(z{i,1})).m(L.MT(2)).xs(:,:,L.T(t),L.S0(b)) = ...
                    lib.elScatKernel(z{i,2},t);
                L.z(L.ZAID(z{i,1})).m(L.MT(251)).xs(:,L.T(t),L.S0(b)) = ...
                    lib.muScat(z{i,2},t);
            end
        end
    end

elseif ~isstruct(p.makeLibraryTempFlag) && p.makeLibraryTempFlag == 2
    % This is the preferred mode of operation currently (not long term).

    % Load previously-generated cross sections from Geoff's or Andy's MAT
    % files.
    load('FULLLIBRARY.mat');
    load('water.mat');
    load('reducedlib.mat');
    load('NAxslib.mat');

    % Subfunction below.
    L = initLibrary(p);

    % Cell arrays for indexing between old and new cross section data.
    z = {1001, 'H1';
         8016, 'O16';
         92235, 'U235';
         92238, 'U238';
         11023, 'NA23'};

    m = {18,'18';
         102,'102';
         251,'251';
         452,'452'};

    % For each ZAID.
    for i = 1:5
        if z{i,1} == 92235 || z{i,1} == 92238
            L.z(L.ZAID(z{i,1})).isFissionable = 1;
        end

        for t = L.Ts
            for b = L.S0s

                % Manage the fact that S0 may be -1.
                if b == -1
                    bstr = '_m1';
                else
                    bstr = ['_' num2str(b)];
                end
                if i == 1 || i == 2
                    bstr = '';
                end

                % For each MT that has an S0 dimension.
                for j = 1:2
                    if (j == 1 && L.z(L.ZAID(z{i,1})).isFissionable == 0)
                        % Don't look at fission for nonfissile ZAIDs.
                    else
                        % Assignment statement.
                        eval(['L.z(L.ZAID(z{i,1})).' ...
                              'm(L.MT(m{j,1})).' ...
                              'xs(:,L.T(t),L.S0(b)) = ' ...
                              z{i,2} '_3_' m{j,2} '_' num2str(t) bstr ';']);
                    end
                end

                for j = 3
                    if z{i,1} == 11023
                        eval(['L.z(L.ZAID(z{i,1})).' ...
                              'm(L.MT(m{j,1})).' ...
                              'xs(:,L.T(t),L.S0(b)) = ' ...
                              z{i,2} '_3_' m{j,2} '_' num2str(t) '_10;']);
                    else
                        eval(['L.z(L.ZAID(z{i,1})).' ...
                              'm(L.MT(m{j,1})).' ...
                              'xs(:,L.T(t),L.S0(b)) = ' ...
                              z{i,2} '_3_' m{j,2} '_' num2str(t) ';']);
                    end
                end

                % For each MT that does not have an S0 dimension.
                for j = 4
                    if (j == 4 && L.z(L.ZAID(z{i,1})).isFissionable == 0)
%                      || (j == 1 && z{i,1} == 1001)
                       % Data unavailable.
                    else
                        eval(['L.z(L.ZAID(z{i,1})).' ...
                              'm(L.MT(m{j,1})).' ...
                              'xs(:,L.T(t),L.S0(b)) = ' ...
                              z{i,2} '_3_' m{j,2} '_' num2str(t) ';']);
                    end
                end

                % Elastic scattering kernel. For i == 1 or 2, use _0_2_ for
                % scattering kernel.
                eval(['L.z(L.ZAID(z{i,1})).' ...
                      'm(L.MT(2)).' ...
                      'xs(:,:,L.T(t),L.S0(b)) = ' ...
                      z{i,2} '_' num2str(6*(i > 2)) '_2_' num2str(t) bstr ';']);
            end
        end
    end

    % Take care of water, which knows no S0 right now.
    for temp = L.Ts
        for s0 = L.S0s
                % Manage the fact that S0 may be -1.
                if s0 == -1
                    bstr = '_m1';
                else
                    bstr = ['_' num2str(s0)];
                end

            % Elastic scattering kernel.
            eval(['L.z(L.ZAID(222)).m(L.MT(2)).xs(:,:,L.T(temp),L.S0(s0)) = '...
                  'H2O_0_2_600' ';']);

            % 1-D data (not kernels).
            for j = [2 3]
            eval(['L.z(L.ZAID(222)).m(L.MT(m{j,1})).xs(:,L.T(temp),L.S0(s0)) = '...
                  'H2O_3_' m{j,2} '_600' bstr ';']);
            end

        end
    end
%}
