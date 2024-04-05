function [Ri] = CreatePi42(L, p, r, Ri)
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
%   20120619    cld2469     refactoring away from "fuel" and "moderator"
%
% TASKLIST
%   Taken from Geoff's version of this function:
%   1- Deal with single transmission cross sections.
%   2- Deal with alternate geometries. This could include writing a function
%   that takes geometry and returns x and c for Sauer's formulation.
%
%   3- The current way this is written is as fuel vs moderator. Is that general
%   enough to account for all types of reactors? Pebble-bed?
%   4- Create warning if a cell with fisisonable isotopes is set to the
%   antipin.

import vbudsii.*

Cell = Ri.Cell;

util.PrintEntering(p, 'CreatePi');
disp(' ===========================  Homogeneous 2');
if r.nCells ~= 4
    error('Homogenize only works for 4 cells right now.');
end

% Unit cell definition (uc stands for unit cell).
uc = r.uc;

% I could imagine, here, having a boolean check for uc.type = 'pincell'.

% Determine if cells are fuel or moderator. This is probably better as a user
% input for cylinder vs anticylinder. To do this I need to have a much more
% flexible way for users to assign cells to physically meaningful spaces.
pinidxs = [];
antipinidxs = [];

if isfield(uc, 'AntipinCellNames') && isfield(uc, 'PinCellNames')
    % Count the number of fuel and moderator cells.
    nAntipinCells = length(uc.AntipinCellNames);
    nPinCells = length(uc.PinCellNames);
    
    % For each cell.
    for cellidx = 1:r.nCells
        for iAntipin = 1:nAntipinCells
            if strcmp(r.cellDef(cellidx).name, uc.AntipinCellNames(iAntipin))
                antipinidxs = [antipinidxs cellidx];
            end
        end
        for iPin = 1:nPinCells
            if strcmp(r.cellDef(cellidx).name, uc.PinCellNames(iPin))
                pinidxs = [pinidxs cellidx];
            end
        end
    end
else    
    % Determine if cells are pin or antipin. This is probably better as a user
    % input for cylinder vs anticylinder. To do this I need to have a much more
    % flexible way for users to assign cells to physically meaningful spaces.
    pinidxs = [];
    antipinidxs = [];
    
    % For each cell.
    for cellidx = 1:r.nCells
        if r.cellDef(cellidx).isFissionable
            % The cell contains fissile ZAIDs, and so is a pin cell.
            pinidxs = [pinidxs cellidx];
        else
            % The cell does not contain fissile ZAIDs, and so is a antipin cell.
            antipinidxs = [antipinidxs cellidx];
        end
    end
    
    % Count the number of fuel and moderator cells.
    nAntipinCells = length(pinidxs);
    nPinCells = length(antipinidxs);
end

% Error-check:
%   -there is at least one antipin cell and one pin cell.
%   -no cell is used twice
%   -all cells are used
if nAntipinCells < 1 || nPinCells < 1
    error('Must have at least one pin and one antipin.');
end

for iCell = 1:r.nCells
    if isempty(find(iCell == antipinidxs)) && isempty(find(iCell == pinidxs))
        error('Cell %s is not assigned to an antipin or pin.', ...
            r.cellDef(iCell).name);
    end
    if length(find(iCell == antipinidxs)) > 1 || ...
            length(find(iCell == pinidxs)) > 1
        error('Cell %s is assigned to more than one pin or antipin.', ...
            r.cellDef(iCell).name);
    end
end

% Make the rest of the code aware of the order of the cells: antipin cells go
% first, in the order they were provided, then come pin cells. Each
% implementation of CreatePi.m will need blah got distracted mid-sentence.
% The user can specify cells in any order, but CreatePi has a certain order.
% This here makes other functions aware of this order, and allows the functions
% to access the proper cell in CreatePi.
pi2cellIdxs = [antipinidxs pinidxs];
cell2piIdxs = zeros(1, r.nCells);
cell2piIdxs(pi2cellIdxs) = 1:r.nCells;

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

mtForXS = 8;
if isfield(p, 'totalfortransport') && p.totalfortransport
    mtForXS = 7;
end

% Alright now we calculate the homogeneity as a function of lambda, mean free
% path.
%
% A is the first pair of fuel-mod
% B is the second pari

if any(uc.f ~= uc.g)
    error('volume ratios must be the same for this method.');
end
%iAPin = 1;
%iAAnti = 2;
%iBPin = 3;
%iBAnti = 4;
%mfpA = 1 ./ Cell(iAAnti).fine(L.MT(mtForXS)).value;
%mfpB = 1 ./ Cell(iBAnti).fine(L.MT(mtForXS)).value;
%% f and g were already checked for normalization
%avgMfp = uc.f(1) * mfpA + uc.f(2) * mfpB;
%% This is important: it gives the center of the erf fuction. Will want to play
%% around with this.
thresh = p.homogenizeThreshFac * mean(uc.pinPitch);
slopeParamAtThresh = p.homogenizeSlopeFac;
%%Homo = 1/4 * (3 + erf(slopeParamAtThresh*(avgMfp - thresh)));
%Homo = 1 * ones(p.nFineGroups, 1);

%figure;
%subplot(2,2,1);
%semilogx(p.fineGroupDef(2:end), [mfpB Cell(iBAnti).fine(L.MT(mtForXS)).value]);
%title('avg mfp vs energy');
%legend('mfpb', 'total mod xs, B');
%subplot(2,2,2);
%semilogx(p.fineGroupDef(2:end), [mfpA mfpB avgMfp]);
%title('should be the same in the case of identical mod.');
%subplot(2,2,3);
%plot(avgMfp, Homo);
%title('homogeneity');
%subplot(2,2,4);
%semilogx(p.fineGroupDef(2:end), [Homo 1 ./ avgMfp])
%title('homogeneity at each energy');
%legend('homogeneity', 'avg mod total xs');
%avgXS = 1 ./ avgMfp;
%save('homogeneityvec', 'Homo', 'avgXS');
%dbstop in vbudsii.multicell.CreatePi4 at 188;

% Antipin cells.
% Sauer's constant is currently defined in the p.uc structure. Information
% about Sauer can be found in the multicell write-ups/documentation held by the
% VBUDSII developers.

if isfield(p, 'useTimsModTransmission')
    for localidx = 1:nAntipinCells
        cellidx = antipinidxs(localidx);
        for iBin = 1:p.nGroups
            fcnHandle = p.useTimsModTransmission;
            TransmitAntipin(localidx, iBin) = fcnHandle( ...
                Cell(cellidx).fine(L.MT(mtForXS)).value(iBin), ...
                uc.pinPitch(localidx), uc.pinDiam(localidx));
        end
        EscapeAntipin(localidx, :) = (1 - TransmitAntipin(localidx,:)) ./ ...
            xAntipin(localidx,:);
        % For the sake of comparison.
        xForCompare = Cell(cellidx).fine(L.MT(mtForXS)).value' * ...
            uc.pinDiam * ...
            (4 * uc.pinPitch(localidx)^2 / (pi * uc.pinDiam(localidx)^2) - 1);
        transmForCompare = 1 ./ ...
            (1 + xForCompare / c ) .^ c;
        Region(regidx).Cell(cellidx).transm = ...
            [TransmitAntipin(localidx, :); ...
            transmForCompare];
    end
else
    c = uc.sauerConst.mod;
    for localidx = 1:nAntipinCells
        % Find the cellidx corresponding to the localidx-th antipinidx so that the
        % Cell structure can be properly indexed.
        cellidx = antipinidxs(localidx);
        xAntipin(localidx, :) = Cell(cellidx).fine(L.MT(mtForXS)).value' * ...
            uc.pinDiam(localidx) * ...
            (4 * uc.pinPitch(localidx)^2 / (pi * uc.pinDiam(localidx)^2) - 1);
        TransmitAntipin(localidx, :) = 1 ./ ...
            (1 + xAntipin(localidx,:) / c ) .^ c;
        %TransmitAntipin(localidx, :) = carlvikrational(xAntipin(localidx,:));
        %disp('WARNING: Modified transmission term.');
        EscapeAntipin(localidx, :) = (1 - TransmitAntipin(localidx,:)) ./ ...
            xAntipin(localidx,:);
    end
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
    cellidx = pinidxs(localidx);
    xPin(localidx, :) = Cell(cellidx).fine(L.MT(mtForXS)).value' * ...
        uc.pinDiam(localidx);
    TransmitPin(localidx, :) = 1 ./ (1 + xPin(localidx, :) / c) .^ c;
    %TransmitPin(localidx, :) = carlvikrational(xPin(localidx,:));
    %disp('WARNING: Modified transmission term.');
    EscapePin(localidx, :) = (1 - TransmitPin(localidx,:)) ./ xPin(localidx,:);
end

% See above discussion for uc.f.
TransmitPinAvg = uc.g * TransmitPin;

% infiniteSum: infinite sum of geometric series to account for
% "telescoping" probability.
denominator = 1 - TransmitAntipinAvg .* TransmitPinAvg;
infiniteSumAntipin = (1 - TransmitAntipinAvg) ./ denominator;
infiniteSumPin = (1 - TransmitPinAvg) ./ denominator;

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
            TransmitPinAvg(binidx) * infiniteSumAntipin(binidx);

    M2 = diag(1 - EscapePin(:,binidx)) + ...
            diag(uc.g) * ones(length(uc.g)) * ...
            diag(EscapePin(:,binidx)) * ...
            TransmitAntipinAvg(binidx) * infiniteSumPin(binidx);

    M3 = diag(uc.g) * ones(length(uc.g), length(uc.f)) *...
            diag(EscapeAntipin(:,binidx)) * infiniteSumPin(binidx);

    M4 = diag(uc.f) * ones(length(uc.f), length(uc.g)) *...
            diag(EscapePin(:,binidx)) * infiniteSumAntipin(binidx);

    PI(:,:,binidx) = [M1, M4 ; M3, M2];
end

% heterogeneous case.
% 1: need to swap around the terms from geoff's case
% do not use tau or any of that stuff.
iAAnti = 1;
iBAnti = 2;
iAPin = 3;
iBPin = 4;
iA = 1;
iB = 2;

mfpAPin = 1./Cell(1).fine(L.MT(7)).value;
mfpAAnti = 1./Cell(2).fine(L.MT(7)).value;
mfpBPin = 1./Cell(3).fine(L.MT(7)).value;
mfpBAnti = 1./Cell(4).fine(L.MT(7)).value;

HomoAAnti = homo(mfpAAnti, thresh, slopeParamAtThresh);
HomoBAnti = homo(mfpBAnti, thresh, slopeParamAtThresh);
HomoAPin = homo(mfpAPin, thresh, slopeParamAtThresh);
HomoBPin = homo(mfpBPin, thresh, slopeParamAtThresh);

Ri.Cell = Cell;
Ri.pi2cellIdxs = pi2cellIdxs;
Ri.cell2piIdxs = cell2piIdxs;
Ri.Homo = [HomoAPin HomoAAnti HomoBPin HomoBAnti];

for iBin = 1:length(p.nFineGroups)
    % AAnti
    PI(2:4,iAAnti,iBin) = PI(2:4,iAAnti,iBin) * HomoAAnti(iBin);
    scaleUp1 = (1 - HomoAAnti(iBin) * sum(PI(2:4,iAAnti,iBin))) / ...
        sum(PI(1:3,iAAnti,iBin));
    PI(1:3,1,iBin) = scaleUp1 * PI(1:3,iAAnti,iBin);

    % BAnti
    PI(1:3,iBAnti,iBin) = PI(1:3,iBAnti,iBin) * HomoBAnti(iBin);
    scaleUp2 = (1 - HomoBAnti(iBin) * sum(PI(1:3,iBAnti,iBin))) / ...
        sum(PI(2:4,iBAnti,iBin));
    PI(2:4,iBAnti,iBin) = scaleUp2 * PI(2:4,iBAnti,iBin);

    % APin
    PI(2:4,iAPin,iBin) = PI(2:4,iAPin,iBin) * HomoAPin(iBin);
    scaleUp3 = (1 - HomoAPin(iBin) * sum(PI(2:4,iAPin,iBin))) / ...
        sum(PI(1:3,iAPin,iBin));
    PI(1:3,iAPin,iBin) = scaleUp3 * PI(1:3,iAPin,iBin);

    % BPin
    PI(1:3,iBPin,iBin) = PI(1:3,iBPin,iBin) * HomoBPin(iBin);
    scaleUp4 = (1 - HomoBPin(iBin) * sum(PI(1:3,iBPin,iBin))) / ...
        sum(PI(2:4,iBPin,iBin));
    PI(2:4,iBPin,iBin) = scaleUp4 * PI(2:4,iBPin,iBin);
end


%save('heypi', 'PI', 'PIhomo', 'PIhetero');

% Verify that all columns of PI add up to 1. This is conveniently done with the
% sum across dimension 1.
if any(sum(PI, 1)) ~= 1
    sum(PI, 1)
    error('In region %s, columns of PI do not sum to 1.', r.name);
end

if isfield(p, 'saveProbabilities')
    disp(['Saving probability data to the current directory:' pwd]);
    save(p.saveProbabilities, 'EscapeAntipin', 'EscapePin', ...
                        'TransmitAntipin', 'TransmitPin', ...
                        'xAntipin', 'xPin', ...
                        'infiniteSumAntipin', 'infiniteSumPin');
end

Ri.PI = PI;

util.PrintExiting(p, 'CreatePi');

end


function Homo = homo(avgMfp, thresh, slopeParamAtThresh);
Homo = 1/2 * (1 + erf(slopeParamAtThresh*(avgMfp - thresh)));
end

