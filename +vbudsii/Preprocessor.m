function [L, p, g, Region, fissionSpectrum] = Preprocessor(p, g, varargin)
%
%VBUDSII/PREPROCESSOR/PREPROCESSOR Performs error checking on user input, and
% initializes/creates data structures Region, L, and the fission neutron
% spectrum.
%
% INPUT
%   p           Parameter structure.
%   g           Geometry structure.
%
% OUTPUT
%   L           Cross section/data library containing data from ENDF/B and
%               processed by NJOY. The structure is initially created in the
%               MakeLibrary.m function. The script that calls VBUDSII may also
%               contain code that analyzes Results, and accordingly it may be
%               useful to have the Library available for that analysis.
%   p           Parameter structure. Should be the same as the input parameter
%               structure.
%   g           Geometry structure. Should not be the same as the input
%               parameter structure, as the Preprocessor module adds some
%               fields to the structure.
%   Results     Structure containing all results, including Region structure,
%               and the program's runtime.
%   fissionSpectrum     [p.nFineGroups x 1 double] Fission neutron spectrum.
%
% DEPENDENCIES
%   VBUDSII/DATAPROCESSING/MAKELIBRARY
%   VBUDSII/DATAPROCESSING/INTERPFISSIONSPECTRUM

% NOTES
%   The user decides if the cross section library should be remade or not. It
%   must be remade if the current library is based on a fine group structure
%   from that which the user wants to use. This module is also the place where
%   inputs are checked for some easy to detect errors. The geometry structure g
%   is modified so that it knows which regions and cells have fissile ZAIDs.
%   This may not be necessary, and is more of a time efficiency thing. If it is
%   implemented, then this check might need to be redone at each time step.
%   Probably the most important task this module does is initialize the runtime
%   data structure Region. Most fields are initialized to zero, while other
%   fields are given inital values as stipulated in the parameter structure p
%   or the geometry structure g.
%
% MAJOR REVISIONS
%   date        handle      description
%   20111030    cld2469     writing comments
%
% TASKLIST
%   1- Complete error checking.
%   2- Removed the flexibility of different group structures (would need NJOY
%   to change this anyway).

import vbudsii.*

util.PrintEntering(p, 'Preprocessor');

% Get the vbudsii root directory.
S = dbstack('-completenames');
vbudsiiroot = fileparts(S(1).file);

% Make cross section library L, based on user input p.makeLibrary.
if p.makeLibrary == 0 & isempty(varargin{1})

    % Load L generated from NJOY scripts.
    load(fullfile(vbudsiiroot, '+data', p.XSLibraryMAT));

    % If user's supplied edges do not match the edges of the library,
    % the library has to be regenerated. EDIT
    %njoyGroupsMatch = sum(p.fineGroupDef ~= njoyGroupStruct) > 0;

    %if ~njoyGroupsMatch
    %
    %   disp(['WARNING: The XSLibrary is being remade from NJOY because '...
    %        'the stored groups and the group definition you provided '...
    %        'do not match.'])
    %
    %    % Overwite L.
    %    L = MakeLibrary(p);
    %
    %end

elseif p.makeLibrary == 1

    L = dataprocessing.MakeLibrary(p);
    % TODO should only save this if the user asks to save it.
    save(fullfile(vbudsiiroot, '+data', p.XSLibraryMAT), 'L', '-v7.3');
elseif ~isempty(varargin{1})
    L = varargin{1};
end

%% ERROR CHECKING

% Error checking for p.
% Verify number of groups.
if p.nFineGroups ~= length(p.fineGroupDef)-1
    error(['The value of p.nFineGroups does not match the number of groups'...
        ' indicated by p.fineGroupDef.']);
end

if p.nFewGroups ~= length(p.fewGroupDef)-1
    error(['The value of p.nFewGroups does not match the number of groups'...
        ' indicated by p.fewGroupDef.']);
end


% Temperatures are integers

% Boolean inputs are actually booleans.

% nTimeSteps is not zero, and should be capped at something.

% Error checking for g.

% Check number of regions and cells.
if g.nRegions ~= length(g.regionDef)
    error(['The value of g.nRegions does not match the number of regions'...
        ' indicated by the number of regions that have been defined.']);
end

for regidx = 1:g.nRegions
    if g.regionDef(regidx).nCells ~= length(g.regionDef(regidx).cellDef)
        error(['The value of g.regionDef(' num2str(regidx) ').nCells does'...
            ' not match the number of cells indicated by the number of'...
            ' cells that have been defined, which is ' ...
            num2str(length(g.regionDef(regidx).cellDef)) '.']);
    end
end

% Check lengths of everything.

% Warning if a region or cell does not have a name.

% These error checks are for multicell/CreatePi.
for regidx = 1:g.nRegions
    uc = g.regionDef(1).uc;
    if abs(1 - sum(uc.g)) >= 0.005
        error(['g values sum to ' num2str(sum(uc.g)) '. Should sum to 1'])
    end
    if abs(1 - sum(uc.f)) >= 0.005
        error(['f values sum to ' num2str(sum(uc.f)) '. Should sum to 1'])
    end
end


%% MODIFYING g

% Populate isFissionable for the cells, from the ZAIDs.
% Check each region.
for regidx = 1:g.nRegions

    % Assume this region contains no fissile isotopes.
    regFissionable = 0;

    % Check each cell.
    for cellidx = 1:g.regionDef(regidx).nCells

        % Assume this cell contains no fissile isotopes.
        cellFissionable = 0;

        % Check each ZAID.
        for zaid = g.regionDef(regidx).cellDef(cellidx).initZAIDs

            if L.z(L.ZAID(zaid)).isFissionable == 1
                % The cell and region containing this ZAID are fissile.
                cellFissionable = 1;
                regFissionable = 1;
            end

        end

        % Assign boolean result to this cell.
        g.regionDef(regidx).cellDef(cellidx).isFissionable = cellFissionable;
    end

    % Assign boolean result to this region.
    g.regionDef(regidx).isFissionable = regFissionable;
end

%% INITIALIZE Region

% Define Region, which contains all runtime data and results, from material
% composition to flux to reaction rates. The details of this structure are not
% given here.
Region = struct('spectralFlux',[],...
                'fewFlux',[],...
                'kInf',[],...
                'kInfall', [], ...
                'relativePower',[],... % scalar
                'PI',[],...
                'powerDensity',[],...
                'few',struct('value',[]),...
                'pi2cellIdxs', [], ...
                'cell2piIdxs', [], ...
                'Homo', [], ...
                'Cell',struct('spectralFlux',[],...
                              'fewFlux',[],...
                              'kInf',0,...
                              'temp',0,...
                              'ZAIDs',[],...
                              'numDensities',[],...
                              'S0s',[],...
                              'SEscapes',[],...
                              'one',struct('RR',[],...
                                           'value',[],...
                                           'z',struct('t',[],...
                                                      's',[], ...
                                                      'RR',[])),...
                              'few',struct('RR', [], ...
                                           'value',[],...
                                           'z',struct('t',[],...
                                                      's',[], ...
                                                      'RR', [])),...
                              'fine',struct('value',[],...
                                            'z',struct('t',[],...
                                                       's',[]))));

% If the parameter structure has a field that defines a constant value for S0,
% use that to initialize the value of S0 for all energies and ZAIDs.
if isfield(p,'constS0')
    constS0 = p.constS0;
else
    constS0 = 10;
end

% Used to initialize the last dimension of the Region().Cell().fine().z().t
% arrays
nTs = length(L.Ts);

% For each region.
for regidx = 1:g.nRegions

    Region(regidx).spectralFlux = ...
        zeros(p.nFineGroups,g.regionDef(regidx).nCells);
    Region(regidx).fewFlux = zeros(p.nFewGroups,1);
    Region(regidx).kInf = 0;
    Region(regidx).relativePower = 0;
    Region(regidx).pi2cellIdxs = zeros(1, g.regionDef(regidx).nCells);
    Region(regidx).cell2piIdxs = zeros(1, g.regionDef(regidx).nCells);

    % For each MT.
    for mt = L.mainMTs

        if mt == 2 || mt == 16 || mt == 6
            % Kernel.
            % Region(regidx).fine(L.MT(mt)).value = zeros(p.nFineGroups);
            Region(regidx).few(L.MT(mt)).value = zeros(p.nFewGroups);

%        elseif ( mt == 18 || mt == 452 || mt == 9 ) && ...
%            (g.regionDef(regidx).isFissionable == 0)
            % Fission related MTs for non-fissile region. Do nothing.

        else
            % All other MTs.
            % Region(regidx).fine(L.MT(mt)).value = zeros(p.nFineGroups,1);
            Region(regidx).few(L.MT(mt)).value = zeros(p.nFewGroups,1);

        end

        Region(regidx).one(L.MT(mt)).value = 0;
        Region(regidx).one(L.MT(mt)).RR = 0;

    end

    % For each cells.
    for cellidx = 1:g.regionDef(regidx).nCells

        Region(regidx).Cell(cellidx).spectralFlux = ...
            zeros(p.nFineGroups, 1);
            %g.regionDef(regidx).cellDef(cellidx).initSpectralFlux;

        Region(regidx).Cell(cellidx).fewFlux = zeros(p.nFewGroups,1);

        Region(regidx).Cell(cellidx).ZAIDs = ...
            g.regionDef(regidx).cellDef(cellidx).initZAIDs;

        Region(regidx).Cell(cellidx).numDensities = ...
            g.regionDef(regidx).cellDef(cellidx).initNumDensities;

        % There is an S0 value for each region, cell, group, and ZAID.
        Region(regidx).Cell(cellidx).S0s = ...
            constS0 * ones(p.nFineGroups,...
                length(g.regionDef(regidx).cellDef(cellidx).initZAIDs));
        Region(regidx).Cell(cellidx).SEscapes = ...
            zeros(p.nFineGroups,...
                length(g.regionDef(regidx).cellDef(cellidx).initZAIDs));

        Region(regidx).Cell(cellidx).kInf = 0;
        Region(regidx).Cell(cellidx).temp = ...
            g.regionDef(regidx).cellDef(cellidx).initTemp;
        Region(regidx).Cell(cellidx).relativePower = 0;

        if isfield(p, 'useTimsModTransmission')
            Region(regidx).Cell(cellidx).transm = ...
                zeros(2, p.nFineGroups);
        end

        % For each MT.
        for mt = L.mainMTs

            if mt == 2 || mt == 16 || mt == 6

                Region(regidx).Cell(cellidx).fine(L.MT(mt)).value = ...
                    zeros(p.nFineGroups);

                Region(regidx).Cell(cellidx).few(L.MT(mt)).value = ...
                    zeros(p.nFewGroups);

%            elseif ( mt == 18 || mt == 452 || mt == 9 ) && ...
%                (g.regionDef(regidx).cellDef(cellidx).isFissionable == 0)
                % Do nothing.
            else

                Region(regidx).Cell(cellidx).fine(L.MT(mt)).value = ...
                    zeros(p.nFineGroups,1);

                Region(regidx).Cell(cellidx).few(L.MT(mt)).value = ...
                    zeros(p.nFewGroups,1);

            end
            Region(regidx).Cell(cellidx).one(L.MT(mt)).value = 0;
        end

        % For each ZAID.
        for zaididx = 1:length(g.regionDef(regidx).cellDef(cellidx).initZAIDs)

            % Obtain the ZAID corresponding to a ZAID index for this cell.
            zaid = g.regionDef(regidx).cellDef(cellidx).initZAIDs(zaididx);

            % For each MT.
            for mt = L.mainMTs

                if mt == 2 || mt == 16 || mt == 6

            Region(regidx).Cell(cellidx).fine(L.MT(mt)).z(zaididx).t ...
                    = zeros(p.nFineGroups,p.nFineGroups,nTs);

            Region(regidx).Cell(cellidx).fine(L.MT(mt)).z(zaididx).s ...
                    = zeros(p.nFineGroups);

            Region(regidx).Cell(cellidx).few(L.MT(mt)).z(zaididx).t ...
                    = zeros(p.nFewGroups,p.nFewGroups,nTs);

            Region(regidx).Cell(cellidx).few(L.MT(mt)).z(zaididx).s ...
                    = zeros(p.nFewGroups);

%                elseif ( mt == 18 || mt == 452 || mt == 9 ) && ...
%                    L.z(L.ZAID(zaid)).isFissionable == 0
                else

            Region(regidx).Cell(cellidx).fine(L.MT(mt)).z(zaididx).t ...
                    = zeros(p.nFineGroups,nTs);

            Region(regidx).Cell(cellidx).fine(L.MT(mt)).z(zaididx).s ...
                    = zeros(p.nFineGroups,1);

            Region(regidx).Cell(cellidx).few(L.MT(mt)).z(zaididx).t ...
                    = zeros(p.nFewGroups,nTs);

            Region(regidx).Cell(cellidx).few(L.MT(mt)).z(zaididx).s ...
                    = zeros(p.nFewGroups,1);
            
            Region(regidx).one(L.MT(mt)).z(zaididx).RR = 0;
            Region(regidx).few(L.MT(mt)).z(zaididx).RR = ...
                zeros(p.nFewGroups, 1);

                end
            Region(regidx).Cell(cellidx).one(L.MT(mt)).z(zaididx).t = ...
                zeros(1,nTs);
            Region(regidx).Cell(cellidx).one(L.MT(mt)).z(zaididx).s = 0;

            end
        end
    end
end

%% Chi
% Chi needs to be interpolaed only once, since it only depends on the fine
% group definition p.fineGroupDef.
fissionSpectrum = dataprocessing.InterpFissionSpectrum(p);

util.PrintExiting(p, 'Preprocessor');

end


%{
%%%% modifying p
% 6 thru 9 are our own MTs for transport, total, nufission, etc.
p.myMTs = [L.mainMTs 6 7 8 9];
p.myMT(p.myMTs) = 1:length(p.myMTs);
p.myMT = sparse(p.myMT);
%%%% modifying p
%}
%{
% Total number of cells. I'm pretty sure this is never needed.
nTotalCells = 0;
for regidx = 1:g.nRegions
    nTotalCells = nTotalCells + g.regionDef(regidx).nCells;
end
%}
