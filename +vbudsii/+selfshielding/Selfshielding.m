function [Ri S0iterConverged] = Selfshielding(L, p, r, Ri)
%
%VBUDSIIA/SELFSHIELDING/SELFSHIELDING Calls SELFSHIELDING/CALCULATES0 for each
% ZAID in each cell, and stores a value of S0 (or updates the value of S0) that
% is stored in Ri for each ZAID in each Cell. Interpolating the cross section
% library for the calculated values of S0 is done in DATAPROCESSING/RESOLVES0
%
% INPUT
%   L           Library structure.
%   p           Parameter structure.
%   r           Field of the geometry structure: g.regionDef(regidx)
%   Ri          Region(regidx): a region in the Region structure.
%
% OUTPUT
%   Ri          Region(regidx): a region in the Region structure.
%   S0iterConverged [boolean] flag indicating if all values of S0 have
%               converged.
%
% DEPENDENCIES
%   VBUDSII/SELFSHIELDING/CALCULATES0

% NOTES
%   S0 does not have an MT dimension (there is only one value of S0 for each
%   ZAID in a cell).
%
% MAJOR REVISIONS
%   date        handle      description
%   20111103    cld2469     writing comments
%
% TASKLIST
%   1- Consider more elaborate convergence criteria.

import vbudsii.*
% Assume that all S0's have converged (S0 for each ZAID in each cell).
S0iterConverged = 1;

% For each cell.
for cellidx = 1:r.nCells

    % For each zaididx.
    for zaididx = 1:length(Ri.Cell(cellidx).ZAIDs)

        % The actual calculation of S0 has been outsourced.
        [S0 xsMicroEscape] = selfshielding.CalculateS0(L, ...
            p, r, Ri.Cell, cellidx, zaididx, ...
            Ri.PI, Ri.cell2piIdxs);

        if abs(Ri.Cell(cellidx).S0s(:,zaididx) - S0)./ S0 > p.S0iterthresh
            % S0 has not converged for this particular ZAID in this cell, based
            % on a convergence threshold specified in the p structure. The Ri
            % structure holds the value of S0 from the previous iteration.
            S0iterConverged = 0;
        end

        % Assign to the structure the new value of S0.
        Ri.Cell(cellidx).S0s(:,zaididx) = S0;
        Ri.Cell(cellidx).SEscapes(:,zaididx) = xsMicroEscape;

    end
end


