function Cell = ResolveTemp(L, p, r, Cell)
%
%VBUDSII/DATAPROCESSING/RESOLVETEMP Interpolates the cross section library L, at
% the input temperature, for all cells, ZAIDs, and MTs in this region.
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
%   The function loops through all cells, ZAIDs, and MTs, and interpolates
%   cross sections according to the reactor's temperature. The temperatures at
%   which cross sections are stored are defined in the NJOY scripts that
%   generate cross sections. Geoff had a bit of code that dealt with
%   funnelling temperature inputs to certain values.
%
% MAJOR REVISIONS
%   date        handle      description
%   20111031    cld2469     writing comments
%   20111106    cld2469     temperature per cell
%
% TASKLIST
%   1- Need to move to a system where temperature is defined per cell, rather
%   than only having one temperature for the entire reactor.

import vbudsii.*

util.PrintEntering(p, 'ResolveTemp');

% For each cell.
for cellidx = 1:r.nCells

    % Grab the temperature setting for this cell.
    tempset = Cell(cellidx).temp;

    % Process the set temperature, adjusting the temperature if necessary. This
    % if statement comes from Geoff's VBUDSII.
    if tempset < L.Ts(1),
        % Cell temperature is less than the minimum temperature at which cross
        % sections are stored in the library.
    
        temp = L.Ts(1);
        disp(['WARNING: Temperature selection is below ' num2str(L.Ts(1))...
            ' K, using ' num2str(L.Ts(1)) ' K instead.']);
    
    elseif tempset > L.Ts(end)
        % Cell temperature is greater than the maximum temperature at which
        % cross sections are stored in the library.
    
        temp = L.Ts(end);
        disp(['WARNING: Temperature selection is above ' num2str(L.Ts(end))...
            ' K, using ' num2str(L.Ts(end)) ' K instead.']);
    
    elseif sum(tempset == L.Ts) == 0 && sum(abs(L.Ts - tempset) <= 5)
        % If the input temperature is not any of the temperates defined in the
        % library and the input temperature is within 5 Kelvin of a
        % temperature defined in the library, use as the temperature the
        % temperature that we are within 5 K of.
    
        temp = L.Ts( find(L.Ts > (tempset-6), 1) );
        disp(['WARNING: Since the reactor temperature ' num2str(temp) ...
            ' is within 5 K of ' num2str(temp) ' K we will use ' ...
            num2str(temp) ' K.']);
    
    else
        temp = tempset;
    end
    
    if any(temp == L.Ts)
        % We hit a temperature spot-on.
        tempidx = L.T(temp);
        % For each zaididx.
        for zaididx = 1:length(Cell(cellidx).ZAIDs)
            
            % Determine ZAID from zaididx.
            z = Cell(cellidx).ZAIDs(zaididx);
            
            % For each MT.
            for m = L.mainMTs
                
                if m == 2 || m == 16 || m == 6 % 16 is n,2n, it's a kernel!!
                    % Scattering kernel.
                    Cell(cellidx).fine(L.MT(m)).z(zaididx).t = ...
                        squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,:,tempidx,:));
                else
                    % Not kernels, etc.
                    Cell(cellidx).fine(L.MT(m)).z(zaididx).t = ...
                        squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,tempidx,:));
                end % m == 2
            end % for m
        end % for zaididx
    else
        % Determine the lower index for temperature (the greatest temperature at
        % which cross sections are stored that is less than the input temperature).
        TLoweridx = max(1, min( find( temp - L.Ts >= 0, 1, 'last'), length(L.Ts)-1));
        % Index in the library of the least temperature that is greater than the
        % input temperature.
        TUpperidx = min(TLoweridx + 1, length(L.Ts));
        
        % Actual temperatures at these indices.
        
        TLower = L.Ts(TLoweridx);
        TUpper = L.Ts(TUpperidx);
        
        % For each zaididx.
        for zaididx = 1:length(Cell(cellidx).ZAIDs)
            
            % Determine ZAID from zaididx.
            z = Cell(cellidx).ZAIDs(zaididx);
            
            % For each MT.
            for m = L.mainMTs
                
                if m == 2 || m == 16 || m == 6 % 16 is n,2n, it's a kernel!!
                    % Scattering kernel.
                    XSL = squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,:,TLoweridx,:));
                    XSU = squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,:,TUpperidx,:));
                    
                    %    elseif ( m == 18 || m == 452 || m == 9 ) && ...
                    %        L.z(L.ZAID(z)).isFissionable == 0
                else
                    
                    % Not kernels, etc.
                    XSL = squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,TLoweridx,:));
                    XSU = squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,TUpperidx,:));
                    
                end % m == 2
                
                % Store interpolation result into Cell structure for output.
                Cell(cellidx).fine(L.MT(m)).z(zaididx).t = ...
                    1/(TUpper-TLower)* ...
                    ( (temp - TLower)*XSU + ...
                    (TUpper - temp)*XSL );
            end % for m
        end % for zaididx
    end
end % for cellidx

util.PrintExiting(p, 'ResolveTemp');

end
