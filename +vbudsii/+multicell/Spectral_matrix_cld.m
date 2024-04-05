function [R, F, Lmat, Rform, Fform, Lform] = ...
    Spectral_matrix_cld(L, p, r, Ri, fissionSpectrum)
%
%VBUDSII/MULTICELL/SPECTRAL_MATRIX_CLD Assembles PI matrix and macroscopic
% cell-level cross sections into matrices to be used in a matrix equation
% solver.
%
% INPUT
%   L           Library structure.
%   p           Parameter structure.
%   r           Field of the geometry structure: g.regionDef(regidx)
%   Ri          Region(regidx): a region in the Region structure.
%   fissionSpectrum     [p.nFineGroups x 1 double] Fission neutron spectrum.
%
% OUTPUT
%   R           [r.nCells * p.nFineGroups square double] Neutron "source" 
%               from scatter.
%   F           [r.nCells * p.nFineGroups square double] Neutron source from
%               fission.
%   Lmat        [r.nCells * p.nFineGroups x 1 double] Neutron sinks, LHS of
%               spectral equation. Renamed as Lmat because of naming conflict
%               with the Library structure.
%   Rform       [r.nCells x p.nFineGroups double] R multiplied by flux
%               (and reshaped).
%   Fform       [r.nCells x p.nFineGroups double] F multiplied by flux
%               (and reshaped).
%   Lform       [r.nCells x p.nFineGroups double] Lmat multiplied by flux
%               (and resshaped).
%
% R, F, and Lmat obey the following neutron conservation equation:
%  L * phi = R * phi + 1/k_{inf} * F * phi
%
% DEPENDENCIES

% NOTES
%   I do not fully understand a fair amount of the matrix operations that occur
%   in this function.
%
% MAJOR REVISIONS
%   date        handle      description
%   20111106    cld2469     writing comments
%
% TASKLIST
%   1- Slowly move to a different function structure so not so many arguments
%   are being passed back and forth.
%   2- It looks like Geoff wanted this to become part of SPECTRAL_SOLV.
%   3- Purpose of K scalar.

import vbudsii.*

util.PrintEntering(p, 'Spectral_matrix_cld');

cellidxs = Ri.pi2cellIdxs;

% This constant is used in the assignment of the Fout matrix further down. I am not sure of its
% purpose.
K = 1;

%% R matrix - Neutron "source" from scatter.
% R comes from a product of cell volumes and the elastic scattering kernel.

% Initialize preliminary R matrix.
R3 = zeros(r.nCells, p.nFineGroups, p.nFineGroups);

% For each cell.
for piidx = 1:r.nCells
    R3(piidx,:,:) = r.relVolumes(cellidxs(piidx)) * ...
        Ri.Cell(cellidxs(piidx)).fine(L.MT(2)).value;
end;

% Rearrange the matrix R3 and the PI matrix.
R2 = size32(R3, [3, 1, 2], 1);
P1 = size32(Ri.PI, [2, 3, 1], 2);

%I1 will make R(i,u)=>R(u,i). Can I use the function reshape or rot here maybe
%transpose? Chris does not understand what is happening here
I1 = zeros(r.nCells * p.nFineGroups);

for j = 1:(r.nCells * p.nFineGroups)
    rr = ((j - mod(j - 1, r.nCells)) - 1) / r.nCells;
    I1(j,(j - 1) * p.nFineGroups + 1 + rr - rr * r.nCells * p.nFineGroups) = 1;
end

% Finally assign to R.
R = I1 * P1 * R2;

%% F matrix - Neutron source from fission.
% F comes from a product of cell volumes, the fission spectrum, group width in
% legarthy, and the fission cross section.

% 2.303 is the lethargy constant, or ln(10)
if isfield(p, 'weightBydE') && ~p.weightBydE 
    du = log(10) * ones(p.nFineGroups, 1);
    binWidth = du;
else
    % TODO
    disp('Weighting by dE instead of dLethargy');
    binWidth = diff(p.fineGroupDef);
end

% For each cell.
for piidx = 1:r.nCells
    if r.cellDef(cellidxs(piidx)).isFissionable
        % This cell has fissile isotopes.

        % For each initial neutron energy.
        for binidx = 1:p.nFineGroups % (u')

            % Neutrons from each incoming energy cause the creation of fission
            % neutrons that end up in various energy groups.
            F3(piidx,:,binidx) = r.relVolumes(cellidxs(piidx)) * ...
                fissionSpectrum(:) * ...
                Ri.Cell(cellidxs(piidx)).fine(L.MT(9)).value(binidx);
        end

    elseif r.cellDef(cellidxs(piidx)).isFissionable == 0
        F3(piidx,:,:) = zeros(p.nFineGroups);
    end % if r.cellDef(cellidx).isFissionable
    
end % for cellidx

F2 = size32(F3, [3,1,2], 1);
F = I1 * P1 * F2;

%% L matrix - Neutron sink from all reactions.
% L comes from a product of cell volumes and the total cross section.
LL = zeros(r.nCells, p.nFineGroups);

if isfield(p,'transportfortotal') && p.transportfortotal
    mtForSink = 8;
else
    mtForSink = 7;
end
% For each cell.
for piidx = 1:r.nCells
    LL(piidx,:) = r.relVolumes(cellidxs(piidx)) * ...
    Ri.Cell(cellidxs(piidx)).fine(L.MT(mtForSink)).value';
    % The cross section vectors have dimensions of [p.nFineGroups x 1].
end

% Reshape.
Lmat = diag(reshape(LL, r.nCells * p.nFineGroups, 1));

% Copy the flux in the reactor from the last time step into a format that is
% formidable for the matrices being defined in this function.
regionFlux = ones(r.nCells, p.nFineGroups);

% For each cell.
for piidx = 1:r.nCells
    regionFlux(piidx,:) = Ri.Cell(cellidxs(piidx)).spectralFlux';
end

% Intermediate step for computing *form matrices: multiply by reshaped flux.
regionFluxReshaped = reshape(regionFlux, (r.nCells * p.nFineGroups), 1);
Rout = R * regionFluxReshaped;
Fout = F * regionFluxReshaped * K;
Lout = Lmat * regionFluxReshaped;

% Assign *form matrices.
Rform = reshape(Rout, r.nCells, p.nFineGroups);
Fform = reshape(Fout, r.nCells, p.nFineGroups);
Lform = reshape(Lout, r.nCells, p.nFineGroups);

util.PrintExiting(p, 'Spectral_matrix_cld');

end


function [out] = size32(input, order, type)
%This program scales a tensor of the form
%E(i,j,k) into a tensor of the form (diagonalized in the 2nd dimension).
%
%           i=1    i=2    i=3   i=4
%       [ 1     | 1     |               ]
%       [  2    |  2    |               ]
% out=  [   3   | j=3   |               ] k=1
%       [    4  |    4  |               ]
%       [     5 |     5 |               ]
%        ----------------------------
%       [ 1     | 1     |               ] k=2
% etc

%A typical input would take the form
%[M_out]=size32( permute(M_in,[3,2,1]));
%
% type needs to be defined properly

%Helpful functions
%    permute     - Permute array dimensions.
%    ipermute    - Inverse permute array dimensions.
%    shiftdim    - Shift dimensions.
%    reshape     - Change size

X0 = permute(input, order);

[i, j, k] = size(X0);
X1 = reshape(X0, i, j * k)';
if type == 1,
    X2 = [ X1(:,1), zeros(j * k, j - 1)];
    for n = 2:i
        X2 = [ X2, X1(:, n), zeros(j * k ,j - 1)];
    end
    X3 = X2;
    % Chris cannot tell if the previous line is just initializing X3 or if X3
    % ends up with some of the data from X2.
    for n = 1:(j * k)
        r = mod(n - 1, j);
        X3(n,:) = circshift(X3(n,:), [0, r]);
    end
    out = X3;
elseif type == 2,
    for n = 1:j * k
        r = mod(n - 1, j);
        X2(n,:) = [zeros(1, i * r), X1(n,:), zeros(1, i * (j - r - 1))];
    end
    out=X2;
    
end

end




%{
for cellidx = 1:r.nCells
    for binidx = 1:p.nFineGroups
%        LL(:,binidx)=r.relVolumes.*Sigma.T(:,binidx);
        LL(cellidx, binidx) = r.relVolumes(cellidx) * ...
        Cell(cellidx).fine(L.MT(7)).value(binidx);
    end
end
%}

%   F31(:,:)=F3(1,:,:);
%   F32(:,:)=F3(2,:,:);
%
%   F31-2*F32

%% (n,2n) or other absorption
% This section is left blank and could be used at a future time to include (n,2n) reactions (or other absorption)
% It is commented out, but could easily be used, argin needs to include
% N_2N(r.nCells)
%
%       du=2.303*ones(p.nFineGroups,1);  'what is du?'
%
%            for j=1:1:r.nCells
%             if N_2N(j)==1,
%                for up=1:1:p.nFineGroups %(u')
%                 A3(j,:,up)=r.relVolumes(j)*du(up)*Sigma.?(j,up);
%                end;
%             elseif N_2N(j)==0,
%                 A3(j,:,:)=zeros(p.nFineGroups,p.nFineGroups);
%             else, 'ERROR bad Region encoding in region j', j
%             end,
%            end;
%         A2=size32(A3,[3,1,2],1);
%         A=I1*P1*F2;
%
%         R=R+A;

%
% Assumed input forms.
%   Sigma.T(i,u)
%   Sigma.S(j,u,u')   (u <- u')
%   Sigma.F(j,u')
%   Chi(u)
%   PI(i,j,u)         (i <- j )
%   r.relVolumes(j)
%
%   j=1:1:r.nCells - spatial
%   u=1:1:p.nFineGroups    - energy
%
%   Regions(j)      - 1 if fissile matrial
%                   - 0 if no fissile material
%
%   phi input can be removed later
%    phi(i,u)
%
