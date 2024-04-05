function PI = CreatePi(L, p, r, Cell)
%
%VBUDSII/MULTICELL/CREATEPI Generates the PI probability matrix that describes
% the probabilities concerning the movement of neutrons from cell to cell. This
% method can manage more than two cells, compared to other methods that can only
% manage two cells.
%
% INPUT
%   L           Cross section/data library containing data from ENDF/B and
%               processed by NJOY. The structure is initially created in the
%               MakeLibrary.m function. The script that calls VBUDSII may also
%               contain code that analyzes Results, and accordingly it may be
%               useful to have the Library available for that analysis.
%   p           Parameter structure.
%   r           Field of the geometry structure: g.regionDef(regidx)
%   Cell        Field of the Region structure: Region(regidx).Cell
%
% OUTPUT
%   PI          [r.nCells x r.nCells x p.nFineGroups double]
%   Cell        Field of the Region structure: Region(regidx).Cell
%
% DEPENDENCIES

% NOTES
%   1- For Sauer's formula, see Stacey pg 330.
%
% MAJOR REVISIONS
%   date        handle      description
%   20111102    cld2469     writing comments
%
% TASKLIST
%   Taken from Geoff's version of this function:
%   1- Deal with single transmission cross sections.
%   2- Deal with alternate geometries. This could include writing a function
%   that takes geometry and returns x and c for Sauer's formulation.
%
%   3- The current way this is written is as fuel vs moderator. Is that general
%   enough to account for all types of reactors? Pebble-bed?

import vbudsii.*

util.PrintEntering(p, 'CreatePi');

% Unit cell definition (uc stands for unit cell).
uc = p.uc;

% I could imagine, here, having a boolean check for uc.type = 'pincell'.

% Determine if cells are fuel or moderator. This is probably better as a user
% input for cylinder vs anticylinder. To do this I need to have a much more
% flexible way for users to assign cells to physically meaningful spaces.
fuelidxs = [];
modidxs = [];
% For each cell.
for cellidx = 1:r.nCells
    if r.cellDef(cellidx).isFissionable
        % The cell contains fissile ZAIDs, and so is a fuel cell.
        fuelidxs = [fuelidxs cellidx];
    else
        % The cell does not contain fissile ZAIDs, and so is a moderator cell.
        modidxs = [modidxs cellidx];
    end
end

% Count the number of fuel and moderator cells.
nAntipinCells = length(modidxs);
nPinCells = length(fuelidxs);

% Initialize.
% x: "equivalent diameter of the region measured in transport mean free paths".
% T: Transmission probability.
% P: Escape probability.
xAntipin = zeros(nAntipinCells, p.nFineGroups);
TransmitAntipin = zeros(nAntipinCells, p.nFineGroups);
EscapeAntipin = zeros(nAntipinCells, p.nFineGroups);

xPin = zeros(nPinCells, p.nFineGroups);
TransmitPin = zeros(nPinCells, p.nFineGroups);
EscapePin = zeros(nPinCells, p.nFineGroups);

% Antipin cells.
% Sauer's constant is currently defined in the p.uc structure. Information
% about Sauer can be found in the multicell write-ups/documentation held by the
% VBUDSII developers.
c = uc.sauerConst.mod;
for localidx = 1:nAntipinCells
    % Find the cellidx corresponding to the localidx-th modidx so that the Cell
    % structure can be properly indexed.
    cellidx = modidxs(localidx);
    xAntipin(localidx, :) = Cell(cellidx).fine(L.MT(8)).value' * uc.pinDiam * ...
                        (4 * uc.pinPitch^2 / (pi * uc.pinDiam^2) - 1);
    TransmitAntipin(localidx, :) = 1 ./ (1 + xAntipin(localidx,:) / c ) .^ c;
    EscapeAntipin(localidx, :) = (1 - TransmitAntipin(localidx,:)) ./ xAntipin(localidx,:);
end

% This line becomes relevant for regions with more than 1 moderator cell. uc.f
% contains weights, summing to 1, that characterize the relative distribution
% of all the moderator cells. For 1 moderator cell, uc.f = 1.
% Make note that this is actual matrix multiplication, not element-wise
% multiplication.
TransmitAntipinAvg = uc.f * TransmitAntipin;

% Pin cells.
c = uc.sauerConst.fuel;
for localidx = 1:nPinCells
    cellidx = fuelidxs(localidx);
    xPin(localidx, :) = Cell(cellidx).fine(L.MT(8)).value' * uc.pinDiam;
    TransmitPin(localidx, :) = 1 ./ (1 + xPin(localidx, :) / c) .^ c;
    EscapePin(localidx, :) = (1 - TransmitPin(localidx,:)) ./ xPin(localidx,:);
end

% See above discussion for uc.f.
TransmitPinAvg = uc.g * TransmitPin;

% collisionProb: collision probabilities.
denominator = 1 - TransmitAntipinAvg .* TransmitPinAvg;
collisionProbAntipin = (1 - TransmitAntipinAvg) ./ denominator;
collisionProbPin = (1 - TransmitPinAvg) ./ denominator;

%% Create Pi matrix terms.

%   [       |          ]   This matrix has the form PI(i,j,u), where the u
%   [  M1   |    M4    ]   is not visible.  It should be noted that M1 & M2
%   [_______|__________]   are square. M3 & M4 represent the transmission
%   [       |          ]   between mod->fuel, and fuel->mod respectively.
%   [  M3   |    M2    ]
%   [       |          ]   The convention on PI is PI( i <- j , u ).  Thus,
%   [       |          ]   sums are to be on i and u & contractions on j.
% For a derivation of the PI matrix, see the documents of the VBUDSII
% developers.

% Initialize.
PI=ones(r.nCells, r.nCells, p.nFineGroups);

% For each neutron energy group (or "bin").
for binidx=1:1:p.nFineGroups

    M1 = diag(1 - EscapeAntipin(:,binidx)) + ...
            diag(uc.f) * ones(length(uc.f)) * ...
            diag(EscapeAntipin(:,binidx))*...
            TransmitPinAvg(binidx) * collisionProbAntipin(binidx);

    M2 = diag(1 - EscapePin(:,binidx)) + ...
            diag(uc.g) * ones(length(uc.g)) * ...
            diag(EscapePin(:,binidx)) * ...
            TransmitAntipinAvg(binidx) * collisionProbPin(binidx);

    M3 = diag(uc.g) * ones(length(uc.g), length(uc.f)) *...
            diag(EscapeAntipin(:,binidx)) * collisionProbPin(binidx);

    M4 = diag(uc.f) * ones(length(uc.f), length(uc.g)) *...
            diag(EscapePin(:,binidx)) * collisionProbAntipin(binidx);

    PI(:,:,binidx) = [M1, M4 ; M3, M2];
end

util.PrintExiting(p, 'CreatePi');

end


%% File Description
% This file will make use of the Geometry, and Transmission crosssections
% to formulate the Pi Matricies.  At this point we are using Sauer's
% formula.
%
%  Assumes
%     Sigma.R(i,u)
%
%   Geometries needed (Assume pin geometry )
%       p - pin pitch
%       d - pin diameter
%       f -
%       g -
%
%   Transmission cross sections for each region
% T transmission
% P escape
% tau collision probability
% want to switch from variable to descriptive names.
