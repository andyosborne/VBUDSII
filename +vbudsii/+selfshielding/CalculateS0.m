function [s0 xsMicroEscape] = CalculateS0(L, p, r, Cell, ...
        cellidx, zaididx, PI, cell2piIdxs)

% The probability of leaving this cell. Should also be 1 -
% PI(cellidx,cellidx,:)
% Get rid of "2", it doesn't mean anything: outgoing?
%PIsum = zeros(p.nFineGroups,1);
%for cellidx2 = 1:r.nCells
%    if cellidx2 ~= cellidx
%        % MAKE SURE THAT GEOFF'S MATRIX IS CELLIDX->CELLIDX2
%        PIsum = PIsum + squeeze( PI(cellidx,cellidx2,:) );
%    end
%end
% Probability of leaving this cell.
piIdx = cell2piIdxs(cellidx);
probLeaving = ones(p.nFineGroups, 1) - squeeze(PI(piIdx, piIdx, :));

if isfield(p, 'transportfortotal') && p.transportfortotal
    mtToUse = 8;
else
    mtToUse = 7;
end

% This calculation doesn't really need to be done here. It can be done one
% level up. TODO
xsMacroEscape = probLeaving ./ ...
    ( 1 - probLeaving ) .* Cell(cellidx).fine(L.MT(mtToUse)).value; % Andy - looks like Cell(cellidx).fine(L.MT(mtToUse)).value gives total macro xs for given cell.

xsMicroEscape = xsMacroEscape / Cell(cellidx).numDensities(zaididx);

% heterogeneous part
es0 = xsMacroEscape; %Cell(cellidx).numDensities(zaididx);

% homogeneous part
% do I want to use cross sections from Cell or from L here?
% s0 = s0 + Cell(cellidx).fine(L.MT(7)).value -
% L.z(L.ZAID(zaid)).m(L.MT(7)).xs * Cell(cellidx).numDensities(zaididx)
for zaididx2 = 1:length(Cell(cellidx).ZAIDs)
    if zaididx2 ~= zaididx
        es0 = es0 + Cell(cellidx).fine(L.MT(mtToUse)).z(zaididx2).s ... % s is "sigma", declared in Preprocessor(), set in dataprocessing.CollapseCell(). I.e. here it's the total sigma for the given zaid etc.
            * Cell(cellidx).numDensities(zaididx2);
    end
end

es0 = es0 / Cell(cellidx).numDensities(zaididx);

s0 = log10(es0);
