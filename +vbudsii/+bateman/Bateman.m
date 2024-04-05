function Ri = Bateman(L, p, r, Ri)
% 
% dN/dt = (A\phi + B)N
% where A represents interactions and B represents natural decay.

% Develop N (numDensities).
% Assumes one fuel cell in this region.
% What this means for other parts of the code:
% -The user specifies which ZAIDs are to be tracked.

% What this needs to be able to do.
% -Discard the creation of isotopes that the user does not want to track.
% -Protons array needs to consider isomeric states.
% -Allow for multiple fuel cells (I think we need to process each cell
% separately).

import vbudsii.*;

util.PrintEntering(p, 'Bateman');

global nIsotopes;

% Time interval for integrating / time stepping.
deltaTime = 3600 * 24 * 7 * 2; % two weeks, in seconds.

% Reorganize data from Region structure.
% TODO temporary: index of fuel, iFuel is temporary right now.
iFuel = 2;
% Grab one-group flux for this cell.
flux = Ri.Cell(iFuel).oneFlux; % Doesn't exist yet.
Composition = Ri.Cell(iFuel).numDensities;
ZAIDs = Ri.Cell(iFuel).ZAIDs;

% To manage the reactions we need to have the number of neutrons and protons in
% each nuclide in this cell.
nIsotopes = length(Ri.Cell(iFuel).ZAIDs);
Protons = floor(ZAIDs / 1000)
Nucleons = ZAIDs - Protons * 1000;

% This structure contains the information that's actually useful 

nFissionProducts = 1;

% Initialize the interactions matrix
InteractionsMatrix = zeros(nIsotopes + nFissionProducts);
% Account for the capture interaction, etc...
InteractionsMatrix = Capture(InteractionsMatrix, Ri.Cell(iFuel));
InteractionsMatrix = Fission(InteractionsMatrix, Ri.Cell(iFuel));
DecayMatrix = Gamma(InteractionsMatrix, Ri.Cell(iFuel));
DecayMatrix = Alpha(InteractionsMatrix, Ri.Cell(iFuel));
% ... and so forth. Make sure this is extendable / modular.

% Combine the interactions and decay into one matrix.
bigMatrix = flux * InteractionsMatrix + DecayMatrix;

% Integrate using one of the available methods below (in subfunctions).
CompositionNext = Step(bigMatrix, Composition, deltaTime);

Ri.Cell(iFuel).numDensities = CompositionNext;
clear nIsotopes;

end

%% Subfunctions

% Utility.
function nucleonNumber = nucleonNumber(ZAID)
    % TODO this is not efficient, since the floor has to be caled for the
    % protons.
nucleonNumber = ZAIDs - floor(ZAID / 1000) * 1000;
end

function protonNumber = protonNumber(ZAID)
protonNumber = floor(ZAIDs / 1000);
end

function ZAID = zaidIdentifier(protonNumber, nucleonNumber)
ZAID = protonNumber * 1000 + nucleonNumber;
end

function index = indexOfZAID(ZAIDs, ZAID);
% TODO checks for uniquness.
index = find(ZAIDs == ZAID);
if length(index) ~= 1
    error('Duplicate ZAIDs.');
end
end

function truefalse = isotopeNotTracked(ZAIDs, ZAID)
% Returns true if the ZAID is not to be tracked.
% This might need to be more intelligent than just finding ZAID in ZAIDs. e.g.
% fission products.
warning('IsotopeNotTracked is not implemented yet.');
truefalse = find(ZAID == ZAIDs);
end

% Time-stepping (integrating) subfunctions.
function CompositionNext = StepMatrixExp(bigMatrix, Composition, deltaTime)
% This is one method for stepping forward the number densities: matrix
% exponential.

% CompositionNext = expm(ST*bigMatrix*S*Composition*deltaTime);

end

% Interaction subfunctions
function matrix = Capture(matrix, Celli)
% Removes the isotope, and creates an isotope with +1 to the number of
% neutrons.

global nIsotopes;

captureXS = Celli.one(L.MT(102)).z;
% The previous line  might need to be the following.
captureXS = zeros(nIsotopes, 1);
for iZAID = 1:nIsotopes
    captureXS(iZAID) = Celli.one(L.MT(102)).z(iZAID).s;
end

for iZAID = 1:nIsotopes
    inZAID = ZAIDs(iZAID);

    inProtons = protonNumber(thisZAID);
    inNucleons = neutronNumber(thisZAID);

    % This is what makes this the capture interaction.
    outNucleons = inNucleons + 1;
    outZAID = zaidIdentifier(outProtons, outNucleons);
    if isotopeNotTracked(ZAIDs, ZAID)
        error('Isotope not tracked.');
    end

    iIn = indexOfZAID(ZAIDs, inZAID);
    iOut = indexOfZAID(ZAIDs, outZAID);
    matrix(iIn, iIn) = matrix(iIn, iIn) - captureXS(iZAID);
    matrix(iOut, iIn) = matrix(iOut, iIn) + captureXS(iZAID);
end

end

num



util.PrintExiting(p, 'Bateman');

end

%Region = struct('spectralFlux',[],...
%                'fewFlux',[],...
%                'kInf',[],...
%                'relativePower',[],... % scalar
%                'PI',[],...
%                'few',struct('value',[]),...
%                'Cell',struct('spectralFlux',[],...
%                              'fewFlux',[],...
%                              'kInf',0,...
%                              'temp',0,...
%                              'powerDensity',[],...
%                              'ZAIDs',[],...
%                              'numDensities',[],...
%                              'S0s',[],...
%                              'one',struct('RR',[],...
%                                           'value',[],...
%                                           'z',struct('t',[],...
%                                                      's',[])),...
%                              'few',struct('value',[],...
%                                           'z',struct('t',[],...
%                                                      's',[])),...
%                              'fine',struct('value',[],...
%                                            'RR',[],...
%                                            'z',struct('t',[],...
%                                                       's',[]))));
%

%numDensities = zeros(nIsotopes, 1);
%
%for cellidx = 1:r.nCells
%    for mt = L.mainMTs
%        if mt == 2 || mt == 9
%        elseif ((mt == 18 || mt == 452 || mt == 9 ) && ~r.cellDef(cellidx).isFissionable)
%        else
%            Ri.Cell(cellidx).fine(L.MT(mt)).RR = ...
%                sum(Ri.Cell(cellidx).fine(L.MT(mt)).value .*...
%                    Ri.Cell(cellidx).spectralFlux);
%            
%        end
%    end
%end
%
%function isotopeMap = CreateIsotopeMap(Protons, Neutrons)
%
%assert(length(Protons) == length(Neutrons))
%for iZAID = 1:length(Protons)
%    isotopeMap(Neutrons(iZAID), Protons(iZAID)) = iZAID;
%end
%end
