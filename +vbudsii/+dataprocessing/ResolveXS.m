function Ri = ResolveXS(L, p, r, Ri)
%
%VBUDSII/DATAPROCESSING/RESOLVEXS Interpolates cross sections for temperature
% and S0, creates macroscopic cross sections for each cell, creates the PI
% matrix for this region. Then, if set to do so, the module iteratively
% converges on a value of S0 for each cell, each ZAID, each group.
%
% INPUT
%   L           Library structure.
%   p           Parameter structure.
%   g           Geometry structure.
%   Ri          Region(regidx): a region in the Region structure.
%
% OUTPUT
%   Ri          Region(regidx): a region in the Region structure.
%
% DEPENDENCIES
%   VBUDSII/DATAPROCESSING/RESOLVETEMP
%   VBUDSII/DATAPROCESSING/RESOLVES0
%   VBUDSII/DATAPROCESSING/COLLAPSECELL
%   VBUDSII/MULTICELL/CREATEPI
%   VBUDSII/SELFSHIELDING/SELFSHIELDING

% NOTES
%
% MAJOR REVISIONS
%   date        handle      description
%   20111030    cld2469     writing comments
%
% TASKLIST
%   1- Determine the proper way to make cell-wise nu and mubar data.
%   2- Should not need to interpolate by temperature again. Should even be able
%   to just comment out this line inside the loop.
%   3- Perhaps we want a more elaborate S0 convergence criterion.

import vbudsii.*

util.PrintEntering(p, 'ResolveXS');

% Interpolate Library for temperature and S0.
Ri.Cell = dataprocessing.ResolveTemp(L, p, r, Ri.Cell);
Ri.Cell = dataprocessing.ResolveS0(L, p, r, Ri.Cell);

% Calculate macroscopic cross sections.
Ri.Cell = dataprocessing.CollapseCell(L, p, r, Ri.Cell);

% Create the PI matrix.
if isfield(p, 'homogenize') && p.homogenize
%        [Ri.PI Ri.pi2cellIdxs Ri.cell2piIdxs] = ...
        if p.homogenize == 1
            Ri = ...
            multicell.CreatePi4(L, p, r, Ri);
        elseif p.homogenize == 2
            Ri = ...
            multicell.CreatePi42(L, p, r, Ri);
        else
            error('!!!')
        end
else
    [Ri.PI Ri.pi2cellIdxs Ri.cell2piIdxs] = ...
        multicell.CreatePi(L, p, r, Ri.Cell);
end

% If user wants to the code to determine appropriate values of S0, iterate.
keepLoopingS0 = p.resolveXS;

itercount = 0;

while keepLoopingS0

    % Calculate S0 for each cell, ZAID, and group.
    [Ri S0iterConverged] = selfshielding.Selfshielding(L, p, r, Ri);

    % Re-interpolate from the library.
    % Technically there should be a way to not need to interpoalte by
    % temperature again.
    % Ri.Cell = ResolveTemp(L, p, r, Ri.Cell);
    Ri.Cell = dataprocessing.ResolveS0(L, p, r, Ri.Cell);

    % Collapse again.
    Ri.Cell = dataprocessing.CollapseCell(L, p, r, Ri.Cell);

    % Re-calculate PI.
if isfield(p, 'homogenize') && p.homogenize
%        [Ri.PI Ri.pi2cellIdxs Ri.cell2piIdxs] = ...
        if p.homogenize == 1
            Ri = ...
            multicell.CreatePi4(L, p, r, Ri);
        elseif p.homogenize == 2
            Ri = ...
            multicell.CreatePi42(L, p, r, Ri);
        else
            error('!!!')
        end
else
    [Ri.PI Ri.pi2cellIdxs Ri.cell2piIdxs] = ...
        multicell.CreatePi(L, p, r, Ri.Cell);
end

    % Keep looping if S0 has not converged AND we are to resolve XSs for S0.
    keepLoopingS0 = ~S0iterConverged;

    % Perhaps we want an elaborate convergence criterion.
    % S0Converged = SelfshieldingConvergence(Ri.Cell);

    % Update the iteration counter.
    itercount = itercount + 1;

    % Print to the command window the current S0 iteration. Typically there are
    % only two iterations required for convergence.
    if p.verbose
        fprintf('S0 iter count: %d\n',itercount);
    end

end

% Write over certain cross sections so that we're usng MCNPX cross sections.
if isfield(p,'useMCNPXTallyXS')
    fprintf('Using MCNPX cross sections. Rerunning CreatePi.\n');
    TallyP = p.useMCNPXTallyXS;

    % Use, to the extent possible, cross sections taken from MCNPX.
    % This is only temporary, and works for only a certain reactor setup.
    %for cellidx = 1:r.nCells
    %    mcnpxflux = Tally{3-cellidx}.value{1}(1:end-1);
    %
    %    % Replacing elastic scattering, nu-fission, and total.
    %    % Must make scatterig diagonal because VBUDSII holds a matrix/kernel
    %    % for elastic scattering, while MCNPX holds a vector.
    %   % Ri.Cell(cellidx).fine(L.MT(2)).value = ...
    %   %     diag(Tally{5-cellidx}.value{2}(1:end-1) ./ mcnpxflux);
    %    Ri.Cell(cellidx).fine(L.MT(7)).value = ...
    %        Tally{5-cellidx}.value{1}(1:end-1) ./ mcnpxflux;
    %    Ri.Cell(cellidx).fine(L.MT(9)).value = ...
    %        Tally{5-cellidx}.value{5}(1:end-1) ./ mcnpxflux;

    %end

    for cellidx = 1:r.nCells
        % elastic kernel 2: just this for now. SKIP FOR WATER?
        for zidx = 1:length(Ri.Cell(cellidx).ZAIDs)
            zaid = Ri.Cell(cellidx).ZAIDs(zidx);
            if zaid == 222
                % Must do a patch, since in MCNPX water is defined with two
                % ZAIDs, so that zidx = 1 is either 8016 or 1001 as opposed to
                % "microscopic" water.
                if length(Ri.Cell(cellidx).ZAIDs) ~= 1
                    error(['VBUDSII cannot manage a cell that contains ' ...
                        'ZAID 222 along with other ZAIDs']);
                end
                % Go from macroscopic to microscopic cross sections.
                prefactor = 1 / Ri.Cell(cellidx).numDensities(zidx);
                mzidx = 0;
            else
                % For all ZAIDs other than water.
                prefactor = 1;
                mzidx = zidx;
            end

            % Input file "correction."
            if isfield(p, 'iFuelPin')
                iPin = p.iFuelPin;
            else
                iPin = 2;
            end
            if iPin == 2
                iMCNP = 3 - cellidx;
            else
                iMCNP = cellidx;
            end

            if isfield(p, 'addinInelastic') && p.addinInelastic
                for MT = L.MTs(find(L.MTs >= 51 & L.MTs <= 91))
                    TallyP(tallystring(num2str(iMCNP*10+4), ...
                        num2str(iMCNP*100+mzidx), '2')) = ...
                        TallyP(tallystring(num2str(iMCNP*10+4), ...
                        num2str(iMCNP*100+mzidx), '2')) + ...
                        TallyP(tallystring(num2str(iMCNP*10+4), ...
                        num2str(iMCNP*100+mzidx), num2str(MT)));
                end
            end


            for gidx = 1:L.nGroups
% TODO add back our own inelastic.
                % scattering 2

                % the following commented out lines don't consider that
                % MCNPX only gives us a vector and not a kernel.
                %Ri.Cell(cellidx).fine(L.MT(mts{2,1})).value(:,gidx) =
                %TallyP(tallystring(cs{cellidx,2}, '2'));
                % How is this one going to work for water? Just skip?
                thisMCNPXxsvec = prefactor * ...
                    TallyP(tallystring(num2str(iMCNP * 10 + 4), ...
                    num2str(iMCNP * 100 + mzidx), '2'));
                MbyV = thisMCNPXxsvec(gidx) / ...
                    sum(Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s(:,gidx));
                Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s(:,gidx) = ...
                    MbyV * Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s(:,gidx);
            end

            % fission 18
            Ri.Cell(cellidx).fine(L.MT(18)).z(zidx).s = ...
                prefactor * ...
                TallyP(tallystring( ...
                num2str(iMCNP * 10 + 4), ...
                num2str(iMCNP * 100 + mzidx), '18'));

            % nu-fission 9
            Ri.Cell(cellidx).fine(L.MT(9)).z(zidx).s = ...
                prefactor * ...
                TallyP(tallystring( ...
                num2str(iMCNP * 10 + 4), ...
                num2str(iMCNP * 100 + mzidx), '-6.-7'));

            % capture 102
            Ri.Cell(cellidx).fine(L.MT(102)).z(zidx).s = ...
                prefactor * ...
                TallyP(tallystring( ...
                num2str(iMCNP * 10 + 4), ...
                num2str(iMCNP * 100 + mzidx), '102'));

            % total: depends on fission, capture, and scatter
            Ri.Cell(cellidx).fine(L.MT(7)).z(zidx).s = ...
                Ri.Cell(cellidx).fine(L.MT(18)).z(zidx).s + ...
                Ri.Cell(cellidx).fine(L.MT(102)).z(zidx).s + ...
                sum(Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s)';

            % TRANSPORT! depends on 
            Ri.Cell(cellidx).fine(L.MT(8)).z(zidx).s = ...
                Ri.Cell(cellidx).fine(L.MT(7)).z(zidx).s - ...
                Ri.Cell(cellidx).fine(L.MT(251)).z(zidx).s .* ...
                sum(Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s)';
        end
    end
end

if isfield(p,'doctorKernel')
    disp('DOCTORING KERNEL <-----------');
    for cellidx = 1:r.nCells
            % elastic kernel 2: just this for now. SKIP FOR WATER?
        for zidx = 1:length(Ri.Cell(cellidx).ZAIDs)
            % Get original groupDef(highEnergyIndex) value.
            energyIndexIn = p.doctorKernel(2);
            energyIndexOut = p.doctorKernel(1);
            initialValue0 = sum(Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s(:,:,1,1))';
            initialValue1 = initialValue0(energyIndexIn);
    
            Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s = ...
                zeros(110);
            Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s(energyIndexOut, ...
                energyIndexIn) = initialValue1;
        end
    end
end

% 1: check to make sure that the prefactor is necessary.
% 2: set up the plot that i want to make
% 3: add in code for the mcnpxkernel
if isfield(p, 'MCNPXhomewaterkernel')
    if isfield(p, 'doctorKernel')
        error('Cannot use MCNPXhomewaterkernel AND doctorKernel options');
    end
    for cellidx = 1:r.nCells
        % elastic kernel 2: just this for now. SKIP FOR WATER?
        for zidx = 1:length(Ri.Cell(cellidx).ZAIDs)
            zaid = Ri.Cell(cellidx).ZAIDs(zidx);
            if zaid == 222
                Ri.Cell(cellidx).fine(L.MT(2)).z(zidx).s = ...
                    p.MCNPXhomewaterkernel;
            end
        end
    end
end

if (isfield(p,'doctorKernel')) || ...
        isfield(p,'useMCNPXTallyXS') || ...
        isfield(p, 'MCNPXhomewaterkernel')

    % remake the other cross sections, though now it's not happening in the
    % same way as before! ideally we'd get microscopic cross sections from
    % MCNPX.

    % we use 2 7 8 9
    % 7 and 8 depend on 2
    Ri.Cell = dataprocessing.CollapseCell(L, p, r, Ri.Cell);
    % Create the PI matrix again, to properly incorporate our changes.
if isfield(p, 'homogenize') && p.homogenize
%        [Ri.PI Ri.pi2cellIdxs Ri.cell2piIdxs] = ...
        if p.homogenize == 1
            Ri = ...
            multicell.CreatePi4(L, p, r, Ri);
        elseif p.homogenize == 2
            Ri = ...
            multicell.CreatePi42(L, p, r, Ri);
        else
            error('!!!')
        end
else
    [Ri.PI Ri.pi2cellIdxs Ri.cell2piIdxs] = ...
        multicell.CreatePi(L, p, r, Ri.Cell);
end
end

if isfield(p, 'runtimeXSplots') && p.runtimeXSplots
for iCell = 1:r.nCells
    iPlot = 0;
    figure;
    for MT = [7 8 9]
        iPlot = iPlot + 1;
        subplot(3,1,iPlot);
        semilogy(Cell(iCell).fine(L.MT(MT)).value);
    end
end
pause
end

util.PrintExiting(p, 'ResolveXS');

end


