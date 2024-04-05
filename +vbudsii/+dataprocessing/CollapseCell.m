function Cell = CollapseCell(L, p, r, Cell)
%
%VBUDSII/DATAPROCESSING/COLLAPSECELL Creates macroscopic cross sections and
% combines mubar (MT = 251) and nu (MT = 452) for each cell in the current region.
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
%   Cell        Field of the Region structure: Region(regidx).Cell
%
% DEPENDENCIES

% NOTES
%   In Geoff's version of this function, which was originally
%   Spectral_matrix_elements, ths function played a more general role and was
%   part of Multicell. See below tasklist for further comments.
%
% MAJOR REVISIONS
%   date        handle      description
%   20111102    cld2469     writing comments
%
% TASKLIST
%   1- Determine the proper way to make cell-wise nu and mubar data. SOLVED:
%   pre-process nufission and transfer.
%   2- As a leftover from Geoff's VBUDSII, cld2469 has left a comment to make
%   an (n,2n) reaction consideration here.
%   3- Should we be running the "for mt = L.mainMTs" loop for all MTs? Is it
%   relevant for cross sections like inelastic scattering or MT = 452? No, it
%   is not. Need to remove 452 and 251 from being cell-ified.

import vbudsii.*

util.PrintEntering(p, 'CollapseCell');

if isfield(p, 'transscale') && p.transscale && ...
        ~p.immutableMyMTs == 2
    error('need immut 2');
end
if isfield(p, 'transscalewateronly') && p.transsalewateronly && ...
        ~p.immutableMyMTs == 2
    error('need immut 2');
end
if isfield(p, 'macromultNuFission') && p.macromultNuFission && ...
        ~p.immutableMyMTs == 2
    error('need immut 2');
end
if isfield(p, 'usen2nn') && p.usen2nn && ~p.immutableMyMTs
    error('To use usen2nn need immutableMyMTs == 2');
end
if isfield(p, 'stammlertransport') && p.stammlertransport && ~p.immutableMyMTs
    error('To use stammlertransport need immutableMyMTs == 2');
end
if isfield(p, 'nomubarforo16') && p.nomubarforo16...
        ~p.immutableMyMTs == 2
    error('need immut 2');
end
if isfield(p, 'hboundmu') && p.hboundmu && ...
        ~p.immutableMyMTs == 2
    error('need immut 2');
end
if isfield(p, 'macromultTransport') && p.macromultTransport && ...
        ~p.immutableMyMTs == 2
    error('need immut 2');
end

% For each cell.
for cellidx = 1:r.nCells
    
    % Abbreviate for brevity. Perhaps not wise for memory usage.
    C = Cell(cellidx);
    
    % For each zaididx.
    for zaididx = 1:length(C.ZAIDs)
        zaid = C.ZAIDs(zaididx);
        
        % This boolean decides if the following 4 MT's are assembled here (after
        % macroscopic cross sections have been made), or are created in
        % MakeLibrary.
        if p.immutableMyMTs == 2
            %        disp('COMPUND XS IN COLLAPSE');
            
            % Absorption, MT = 6.
            if r.cellDef(cellidx).isFissionable
                absorptionxs = C.fine(L.MT(18)).z(zaididx).s + ... %Jess Comment. Does not consider (n,alpha). This is a problem for B-10
                    C.fine(L.MT(102)).z(zaididx).s;
            else
                absorptionxs = C.fine(L.MT(102)).z(zaididx).s;
            end
            absorptionxs=absorptionxs+C.fine(L.MT(107)).z(zaididx).s;% Jess edit, added in (n,alpha) 2.25.2021
            % Must modify mubar before modifying transport (or scattering!).
            if zaid == 1001 && isfield(p, 'hboundmu') && p.hboundmu
                mubar = zeros(length(L.groupDef)-1, 1);
                ethermal = 4e-3; % eV
                eepitherm = 4;
                boundthermal = find( ethermal < L.groupDef, 1, 'first');
                boundepithermal = find( eepitherm >= L.groupDef, 1, 'last');
                mubar(1:boundthermal) = 0.01*ones(boundthermal, 1);
                mubar(boundthermal+1:boundepithermal) = 0.01 + (.67-.01)/ ...
                    (eepitherm - ethermal) * ( ...
                    L.groupDef(boundthermal+1:boundepithermal)' - ethermal);
                mubar(boundepithermal+1:end) = ...
                    0.67*ones(length(L.groupDef)-boundepithermal-1, 1);
                C.fine(L.MT(251)).z(zaididx).s = mubar;
            end
            
            if zaid == 8016 && isfield(p, 'nomubarforo16') && ...
                    p.nomubarforo16
                C.fine(L.MT(251)).z(zaididx).s = zeros(p.nFineGroups, 1);
            end
            
            if zaid == 92238 && isfield(p, 'nomubarforu238') && ...
                    p.nomubarforu238
                C.fine(L.MT(251)).z(zaididx).s = zeros(p.nFineGroups, 1);
            end
            
            if zaid == 92235 && isfield(p, 'nomubarforu235') && ...
                    p.nomubarforu235
                C.fine(L.MT(251)).z(zaididx).s = zeros(p.nFineGroups, 1);
            end
            
            if isfield(p, 'stammlertransport') && p.stammlertransport
                % Must be careful with the possibility of double-modifying the
                % scattering cross section.
                C.fine(L.MT(6)).z(zaididx).s = ...
                    C.fine(L.MT(2)).z(zaididx).s;
                C.fine(L.MT(2)).z(zaididx).s = ...
                    C.fine(L.MT(2)).z(zaididx).s - ...
                    diag(C.fine(L.MT(251)).z(zaididx).s .* ...
                    sum(C.fine(L.MT(2)).z(zaididx).s)');
%                 figure;
%                     plot(C.fine(L.MT(251)).z(zaididx).s .* ...
%                     sum(C.fine(L.MT(2)).z(zaididx).s)');
            end

            if isfield(p, 'usen2nn') && p.usen2nn
                C.fine(L.MT(2)).z(zaididx).s = ...
                    C.fine(L.MT(2)).z(zaididx).s +  ...
                    2 * C.fine(L.MT(16)).z(zaididx).s; 
             
               %     3 * C.fine(L.MT(17)).z(zaididx).s;
                absorptionxs = absorptionxs -  ...
                    2 * sum(C.fine(L.MT(16)).z(zaididx).s)';
            
                %- ...
               %     3 * sum(C.fine(L.MT(17)).z(zaididx).s)'; 
            end
            
            % Total, MT = 7.
            C.fine(L.MT(7)).z(zaididx).s = ...
                sum(C.fine(L.MT(2)).z(zaididx).s)' + ...
                absorptionxs;
            %C.fine(L.MT(4)).z(zaididx).s  + ...
            
            if isfield(p, 'macromultTransport') && p.macromultTransport
            else
                %                disp('multmacroTransport')
                % Transport, MT = 8.
                C.fine(L.MT(8)).z(zaididx).s = C.fine(L.MT(7)).z(zaididx).s - ...
                    C.fine(L.MT(251)).z(zaididx).s.* ...
                    sum(C.fine(L.MT(2)).z(zaididx).s)';
                if isfield(p, 'transscale')
                    if isfield(p, 'transscalewateronly') && ...
                            p.transscalewateronly
                        if strcmp(r.cellDef(cellidx).name, 'coolant')
                            C.fine(L.MT(8)).z(zaididx).s = ...
                                C.fine(L.MT(8)).z(zaididx).s + ...
                                p.transscale * ( ...
                                C.fine(L.MT(7)).z(zaididx).s - ...
                                C.fine(L.MT(8)).z(zaididx).s);
                        end
                    else
                        C.fine(L.MT(8)).z(zaididx).s = ...
                            C.fine(L.MT(8)).z(zaididx).s + ...
                            p.transscale * ( ...
                            C.fine(L.MT(7)).z(zaididx).s - ...
                            C.fine(L.MT(8)).z(zaididx).s);
                    end
                end
                
            end
            if isfield(p, 'stammlertransport') && p.stammlertransport
                C.fine(L.MT(8)).z(zaididx).s = ...
                    NaN*C.fine(L.MT(8)).z(zaididx).s;
            end
            
            if isfield(p, 'macromultNuFission') && p.macromultNuFission
            else
                %            disp('multmacroNuFission')
                % nu-fission, MT = 9.
                if p.makeRealNuFission == 1
                    C.fine(L.MT(9)).z(zaididx).s = ...
                        C.fine(L.MT(452)).z(zaididx).s.*...
                        C.fine(L.MT(18)).z(zaididx).s;
                else
                    % Altered option for comparision to Geoff's VBUDSII.
                    C.fine(L.MT(9)).z(zaididx).s = ...
                        C.fine(L.MT(18)).z(zaididx).s;
                end
            end
        end
        
        % For each MT.
        for mt = [2 6 7 8 9 18 102 251 452]
            
            %            if ( mt == 18 || mt == 452 || mt == 9) && ...
            %                (r.cellDef(cellidx).isFissionable == 0 || ...
            %                 L.z(L.ZAID(zaid)).isFissionable == 0)
            % If either a cell or a zaid within a fissile cell is not
            % fissile, do not add its weighted microscopic cross section
            % contribution to the macroscopic cross section
            
            % Probably can eliminate, from the above boolean, the Cell's
            % isFissionable flag.
            
            %            else
            % The (zaididx~=1) statement causes the value of
            % C.fine(L.MT(mt)).value to not be assigned recursively if the
            % first ZAID is being added. This allows the use of one
            % assignment statement to prevent adding to
            % C.fine(L.MT(mt)).value its value from the last S0 iteration
            % or timestep.
            % This assignment is of the form
            %   \Sigma_{cellidx} = sum(N_{zaididx}  \sigma_{zaididx}).
            C.fine(L.MT(mt)).value = ...
                (zaididx~=1)*C.fine(L.MT(mt)).value + ...
                C.numDensities(zaididx) * ...
                C.fine(L.MT(mt)).z(zaididx).s;
            %            end
        end
        for mt = [251 452]
        end
    end
    
    % This boolean decides if the following 4 MT's are assembled here (after
    % macroscopic cross sections have been made), or are created in
    % MakeLibrary. THIS STUFF IS GARBAGE.
    if p.immutableMyMTs == 0
        %        disp('COMPUND XS IN COLLAPSE');
        
        % Absorption, MT = 6.
        if r.cellDef(cellidx).isFissionable
            absorptionxs = C.fine(L.MT(18)).value + ...
                C.fine(L.MT(102)).value;
        else
            absorptionxs = C.fine(L.MT(102)).value;
        end
        
        % Total, MT = 7.
        C.fine(L.MT(7)).value = sum(C.fine(L.MT(2)).value)' + ...
            absorptionxs;
        %C.fine(L.MT(4)).value  + ...
        
        % Transport, MT = 8.
        C.fine(L.MT(8)).value = C.fine(L.MT(7)).value - ...
            C.fine(L.MT(251)).value.* ...
            sum(C.fine(L.MT(2)).value)';
        
        % nu-fission, MT = 9.
        if p.makeRealNuFission == 1
            C.fine(L.MT(9)).value = C.fine(L.MT(452)).value.*...
                C.fine(L.MT(18)).value;
        else
            % Altered option for comparision to Geoff's VBUDSII.
            C.fine(L.MT(9)).value = C.fine(L.MT(18)).value;
        end
        
    end
    if isfield(p, 'macromultTransport') && p.macromultTransport
        disp('macromultTransport')
        % Transport, MT = 8.
        C.fine(L.MT(8)).value = C.fine(L.MT(7)).value - ...
            C.fine(L.MT(251)).value.* ...
            sum(C.fine(L.MT(2)).value)';
    end
    
    if isfield(p, 'macromultNuFission') && p.macromultNuFission
        disp('macromultNuFission')
        % nu-fission, MT = 9.
        if p.makeRealNuFission == 1
            C.fine(L.MT(9)).value = C.fine(L.MT(452)).value.*...
                C.fine(L.MT(18)).value;
        else
            % Altered option for comparision to Geoff's VBUDSII.
            C.fine(L.MT(9)).value = C.fine(L.MT(18)).value;
        end
    end
    Cell(cellidx) = C;
end

util.PrintExiting(p, 'CollapseCell');

end


