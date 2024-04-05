function Ri = Spectral_Solv(p, r, Ri, L, Rmat, Fmat, Lmat)
%
%VBUDSII/MULTICELL/SPECTRAL_SOLV Computes a matrix solution for an updated
% reactor flux using input matrices.
%
% INPUT
%   p           Parameter structure.
%   r           Field of the geometry structure: g.regionDef(regidx)
%   Ri          Region(regidx): a region in the Region structure.
%   Rmat        [r.nCells * p.nFineGroups square double] Neutron "source"
%               from scatter.
%   Fmat        [r.nCells * p.nFineGroups square double] Neutron source from
%               fission.
%   Lmat        [r.nCells * p.nFineGroups x 1 double] Neutron sinks, LHS of
%               spectral equation.
%
% OUTPUT
%   Ri          Region(regidx): a region in the Region structure, with the
%               spectralFlux field updated.
%
% DEPENDENCIES
%   VBUDSII/MULTCIELL/SPECTRAL_SOLV

% NOTES
%   The method for solution of spectral flux is to cast the problem as a
%   generalized eigenvalue problem
%   A * phi = 1/k_{inf} * B * phi
%   where A = L - R, and B = F
%
% MAJOR REVISIONS
%   date        handle      description
%   20111107    cld2469     writing comments
%
% TASKLIST
%   1- Slowly move to a different function structure so not so many arguments
%   are being passed back and forth.
%   2- Do we really want to be taking the absolute value of our flux solution?
%   3- Make sure we are calculating power correctly.

import vbudsii.*

util.PrintEntering(p, 'Spectral_Solv');

% Call MATLAB's eig routine to obtain eigenvalues and eigenvectors.
[V, D]=eig(Lmat - Rmat, Fmat);

% Find the index of the first finite eigenvalue. This will be the index of
% 1/k_{inf}. In this case, DIAG is undiagonalizing D.
kidx = find( abs( diag(D) ) ~= Inf );

%% ERROR CHECKING
if isempty(kidx)

    error(['ERROR: no spectral flux eigenvalue solution could ' ...
        'be found for region ' r.name '.']);

elseif length(kidx) >= 2
    % There is more than one eigenvalue that is not infinity.

    % Multiple solutions only warrants a warning.
    disp(['WARNING: Spectral flux has multiple solutions in region ' ...
        r.name '. k has the following vals:']);
    disp(sprintf('%.5f, ', 1./D(kidx,kidx)));

    % I am not sure what the next line does. Perhaps finds the vector with the
    % smallest total value that is also all positive.
    [v, b] = min( sum( abs( V(:,kidx) ), 1 ) - abs( sum( V(:,kidx), 1 ) ) );
    kidx = kidx(b);

end

% Used to inspect the eigenvalues to see if there was more than one solution.
% dbstop in Spectral_Solv.m at 72

%% STORE RESULTS

% Take abs to get rid of any imaginary or negative parts?
% The flux is stored in V(:,kidx).
spectralFlux = abs(reshape( V(:,kidx), [r.nCells, p.nFineGroups]))';

% These next lines should be able to compute the reactor's power. However, I do
% not think it is doing that. This should properly be an integral over phi dV.
EnergyPerFission = 200; % MeV/fission e6 * 1.602e-19; % J/fission
EnergyPerFission = 200e6 * 1.602e-19; % J/fission
%powerDensity =  EnergyPerFission * sum(Fmat * abs(V(:,kidx))) / sum(r.relVolumes);
%power =  EnergyPerFission * sum(Fmat * abs(V(:,kidx))) / sum(r.relVolumes);

% 200,000,000 eV/Fission*Fissions/s -> [W]

% Scale the flux output.
%spectralFlux = spectralFlux * p.powerDensity / power;
%spectralFlux = spectralFlux * powerDensity / power;

piidxs = Ri.cell2piIdxs;
Ri.spectralFlux = spectralFlux(:, piidxs);
% For each cell, store the flux for output.
for cellidx = 1:r.nCells
    Ri.relativePower(cellidx) = EnergyPerFission * ...
        sum(Ri.Cell(cellidx).fine(L.MT(9)).value .* ...
        spectralFlux(:,piidxs(cellidx))) * ...
        r.relVolumes(cellidx);
end
%JESSICA EDIT 3.7.19
%Ri.relativePower=100;



Ri.powerDensity = sum(Ri.relativePower) / sum(r.relVolumes);

%JESSICA EDIT 3.7.19
%Ri.powerDensity=1E-01;
Ri.spectralFlux = Ri.spectralFlux / Ri.powerDensity;

for cellidx = 1:r.nCells
    Ri.Cell(cellidx).spectralFlux = Ri.spectralFlux(:,cellidx);
end


% Store k_{inf}.
Ri.kInf = 1/D(kidx,kidx);
Ri.kInfall = 1./diag(D);

util.PrintExiting(p, 'Spectral_Solv');



end
%{

Cell(cellidx).relativePower = ...
    200e3*1.602e-19*sum(F*abs(V(:,kidx)))/sum(r.relVolumes);
% 200,000keV/Fission*Fissions/s -> [kW]

phi = phi * p.powerDensity / Cell(cellidx).relativePower;
%}



