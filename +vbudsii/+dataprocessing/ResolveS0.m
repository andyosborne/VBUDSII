function Cell = ResolveS0(L, p, r, Cell)
%
%VBUDSII/DATAPROCESSING/RESOLVES0 Interpolates the cross section library L, at
% the given values of S0 for all cells, ZAIDs, MTs, and neutron energy groups
% in this region.
%
% INPUT
%   L           Library structure.
%   p           Parameter structure.
%   r           Field of the geometry structure: g.regionDef(regidx)
%   Cell        Field of the Region structure: Region(regidx).Cell
%
% OUTPUT
%   Cell        Field of the Region structure: Region(regidx).Cell
%
% DEPENDENCIES

% NOTES
%
% MAJOR REVISIONS
%   date        handle      description
%   20111031    cld2469     writing comments
%
% TASKLIST
%   1- I am not sure if I am interpolating the scattering kernel correctly: I
%   am currently interpolating with respect to incident energy level.
%   2- Properly initialize XSL and XSU.

import vbudsii.*

util.PrintEntering(p, 'ResolveS0');

% Initialize array of S0 indices into the cross section library.
S0Loweridx = zeros(p.nFineGroups,1);
S0Upperidx = zeros(p.nFineGroups,1);

% Length of the L.S0s list minus 1.
lengthS0min1 = length(L.S0s)-1;

% For each cell.
for cellidx = 1:r.nCells

    % For each zaididx.
    for zaididx = 1:length(Cell(cellidx).ZAIDs)

        z = Cell(cellidx).ZAIDs(zaididx);
        
        % Obtain values of S0 for interpolation.
        S0s = Cell(cellidx).S0s(:,zaididx);

        % Imaginary numbers. This must come before the next error check.
        if ~isreal(S0s)
            disp(['ERROR: S0 has imaginary component for cell ' ...
                   r.cellDef(cellidx).name ' and ZAID ' ...
                   num2str(z) '.']);
            S0s = real(S0s);
        end
        
        % Error-checking if S0 is out of range.
        if p.verbose
            disp(['S0 report cell ' ...
                   r.cellDef(cellidx).name ' ZAID ' ...
                   num2str(z) ' min ' num2str(min(S0s)) ...
                    ' max ' num2str(max(S0s)) '.']);
        end
        if (sum(S0s > max(L.S0s)) > 0) || (sum(S0s < min(L.S0s)) > 0)
            disp(['ERROR: S0 out of range. cell ' ...
                   r.cellDef(cellidx).name ' ZAID ' ...
                   num2str(z) ' min ' num2str(min(S0s)) ...
                    ' max ' num2str(max(S0s)) '.']);
            S0s = min( max( S0s, min(L.S0s)), max(L.S0s));
        end

        if sum(isnan(S0s)) > 0
            disp(sprintf(['ERROR: S0 has NaN values. Cell %s and ' ...
                'ZAID %i.'], r.cellDef(cellidx).name, z));
        end

        % For each neutron energy group, get the greatest S0 for whch S0 is
        % less than the S0 at this cell, ZAID, MT, and group.
        for binidx = 1:p.nFineGroups
            if isnan(S0s(binidx))
                % hack! EDIT.
                S0s(binidx) = L.S0s(2);
            end
        S0Loweridx(binidx) = min(find(S0s(binidx)-L.S0s >= 0, 1, 'last'), ...
            lengthS0min1);

        end
        S0Upperidx = S0Loweridx + 1;
        
        % S0 values at the indicies obtained.
        S0Lower = L.S0s(S0Loweridx)';
        S0Upper = L.S0s(S0Upperidx)';


        % For each MT.
        for m = L.mainMTs
            
            if m == 2 || m == 16 || m == 6 % 16 is n,2n, it's a kernel!!
            % DO NOT KNOW HOW TO HANDLE SCATTERING !!! 
        
                XSL = zeros(p.nFineGroups); 
                XSU = zeros(p.nFineGroups); 

                % For each group.
                for binidx = 1:p.nFineGroups

XSL(:,binidx) = ...
   Cell(cellidx).fine(L.MT(m)).z(zaididx).t(:,binidx,S0Loweridx(binidx));
XSU(:,binidx) = ...
   Cell(cellidx).fine(L.MT(m)).z(zaididx).t(:,binidx,S0Upperidx(binidx));

% Assign result to Cell structure.
Cell(cellidx).fine(L.MT(m)).z(zaididx).s(:,binidx) = ...
    1/(S0Upper(binidx)-S0Lower(binidx)) * ...
    ( (S0s(binidx) - S0Lower(binidx))*XSU(:,binidx) + ...
      (S0Upper(binidx) - S0s(binidx))*XSL(:,binidx));

                end

%            elseif ( m == 18 || m == 452 || m == 9 ) && ...
%                L.z(L.ZAID(z)).isFissionable == 0
            else
                XSL = zeros(p.nFineGroups,1); 
                XSU = zeros(p.nFineGroups,1); 
                % For each group.
                for binidx = 1:p.nFineGroups

XSL(binidx) = ...
    Cell(cellidx).fine(L.MT(m)).z(zaididx).t(binidx,S0Loweridx(binidx));
XSU(binidx) = ...
    Cell(cellidx).fine(L.MT(m)).z(zaididx).t(binidx,S0Upperidx(binidx));

Cell(cellidx).fine(L.MT(m)).z(zaididx).s(binidx) = ...
    1/(S0Upper(binidx)-S0Lower(binidx)) * ...
    ( (S0s(binidx) - S0Lower(binidx))*XSU(binidx) + ...
      (S0Upper(binidx) - S0s(binidx))*XSL(binidx));

                end % for binidx
            end % if m == 2
        end % for m
        clear XSL XSU;
    end % for zaididx
end % for cellidx

util.PrintExiting(p, 'ResolveS0');

end


%XSL = squeeze(Cell(cellidx).fine(g.myMT(m)).z(zaididx).value(:,:,S0Loweridx));
%XSL = squeeze(Cell(cellidx).fine(g.myMT(m)).z(zaididx).value(:,:,S0Upperidx));
