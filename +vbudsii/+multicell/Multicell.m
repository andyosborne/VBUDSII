function Ri = Multicell(L, p, r, Ri, fissionSpectrum)
%
%VBUDSII/MULTICELL/MULTICELL Solves a sources/sinks matrix balance equation
% using the PI matrix.
%
% INPUT
%   L           Library structure.
%   p           Parameter structure.
%   r           Field of the geometry structure: g.regionDef(regidx)
%   Ri          Region(regidx): a region in the Region structure.
%   fissionSpectrum     [p.nFineGroups x 1 double] Fission neutron spectrum.

%
% OUTPUT
%   Ri          Region(regidx): a region in the Region structure, with the
%               spectralFlux field updated.
%
% DEPENDENCIES
%   VBUDSII/MULTICELL/SPECTRAL_MATRIX_CLD
%   VBUDSII/MULTCIELL/SPECTRAL_SOLV

% NOTES
%   Currently the line between the RESOLVEXS and MULTICELL modules is not
%   clearly defined. CREATEPI is a core function of Multicell, but it isn't
%   used here.
%
% MAJOR REVISIONS
%   date        handle      description
%   20111106    cld2469     writing comments
%
% TASKLIST
%   1- Slowly move to a different function structure so not so many arguments
%   are being passed back and forth.

import vbudsii.*

util.PrintEntering(p, 'Multicell');

%if isfield(p,'useMCNPXTallyXS')
    % Use, to the extent possible, some of the cross sections taken from MCNPX.
    % This is only temporary, and works for only a certain reactor setup.
    %Tally = p.useMCNPXTallyXS;
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
    
%end

% Create the matrices involved in the matrix balance equation.
[RM, FM, LM, RMform, FMform, LMform] = ...
    multicell.Spectral_matrix_cld(L, p, r, Ri, fissionSpectrum);

% Solve the matrix equation as an eigenvalue problem.
Ri = multicell.Spectral_Solv(p, r, Ri, L, RM, FM, LM);

util.PrintExiting(p, 'Multicell');

end
