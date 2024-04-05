function [xsOut, fluxOut, groupStructureOut] = CollapseEnergy(groupStructureIn, groupStructureOut, xsIn, fluxIn)

global xstypes;

% EnergyCollapse collapses cross section data from an initial energy group
% structure to a final energy group structure by flux-averaging the cross
% section input. We expect that the input cross section and flux produce
% the same reaction rate as the output cross section and flux.
% EnergyCollapse knows no geometry or material properties. Can process
% multiple cross sections at once if xsIn is formatted as a struct.
%
% Purpose: take 110-group flux output of multicell (for each region of the
% reactor) to create few-group cross sections to be used in a few-group
% diffusion code.
%
% The function has been developed so as to be flexible with application, but in
% its intended use, xsIn is a struct with fields given in cell array
% xstypes, and these input cross sections are microscopic and are for an
% isotope (unit substance). The cross sections come from reducedlib.mat and 
% the flux comes from a multicell code.
%
% INPUT
% groupStructureIn    row vector of energy edges (eV), typically 110 groups
%                     (length 111)
% groupStructureOut   row vector of energy edges (eV), few group. elements must
%                     be a subset of the elements of groupStructureIn
% xsIn                vector or kernel, macro or micro scopic.
%                     If this is a struct, the function performs its action
%                     on fields of the struct given in (global?) cell array
%                     XSTYPES.
% xsIn                size can be vector or kernel, type can be macro or micro
%                     independent of units
% fluxIn              column vector with size of (groupStructureIn-1)
%                     (# cm^-2 s^-1)
%
% All inputs are sorted in increasing energy (not increasing lethargy).
%
% OUTPUT
% xsOut               same datatype, shape and units as xsIn
% fluxOut             column vector with size of (groupStructureOut-1) (# cm^-2 s^-1)
%                     holding rebinned input flux for the corresponding
%                     groupStructureOut, formulated such that
%          sum(xsOut.*fluxOut) = sum(xsIn.*fluxIn)
% groupStructutureOut regurgitation of input groupStructureOut
%
% HOW IT WORKS
%
% Last edited (by and date): cld72@cornell.edu A101018

%% Process inputs

% include global variables (xstypes). not functioning yet
%include_variables

% Deal with some numerical error issue that occurs for the 10^0 = 1 bin.
for iBin = 1:length(groupStructureIn)
    if abs(groupStructureIn(iBin) - 1) < 1e-10
        groupStructureIn(iBin) = 1;
    end
end
% error check, formulated as "this is what's wrong", instead of "this is
% what to fix". Use internal variable names in the messages? what's
% clearer? Talk in the language of the physics/science or of the code?
if length(groupStructureIn) < length(groupStructureOut)
    error('The desired output group structure has more bins than the input group structure.')
end

if sort(groupStructureIn) ~= groupStructureIn
    error('The input group structure is not sorted in ascending order.')
end

if sort(groupStructureOut) ~= groupStructureOut
    error('The output group structure is not sorted in ascending order.')
end

[numRows_GSIn, numCols_GSIn] = size(groupStructureIn);

if numRows_GSIn ~= 1 % "edges" is typically formatted as a row vector.
    error('The input group structure is not in a row vector.')
end

[numRows_GSOut, numCols_GSOut] = size(groupStructureOut);

if numRows_GSOut ~= 1
    error('The output group structure is not in a row vector.')
end

% get the number of input bins and output bins
numBinsIn = numCols_GSIn - 1; % this is also the length (size) of xsIn
numBinsOut = numCols_GSOut - 1; % this is also the length (size) of xsOut

% ensure that groupStructureOut's elements are a subset of groupStructureIn
% elements. first and last element are the same
% other possible checks: no repeated values in either vector, 

% develop the bin indices cell array; will index groupStructureIn
% each cell in GSIidxer is for an output energy bin, and holds the indices
% of the energy bins (not edges) of the input group structure that
% correspond to that bin of the new group structure.
% going from 110 to 3 bins, we have GSIidxer = {1:40,41:90,91:110}
GSIidxer = cell(1,numBinsOut);

for i = 1:numBinsOut % go through each element of the input group structure
     % idx holds the index of the element(s) in groupStructureIn that
     % hold(s) the energy edge groupStructureOut(i);
    isMatchedEdgeL = ( groupStructureIn == groupStructureOut(i) );
    idx1 = find(isMatchedEdgeL); % left edge of this output energy bin
    
    isMatchedEdgeR = ( groupStructureIn == groupStructureOut(i+1) );
    idx2 = find(isMatchedEdgeR); % right edge of this output energy bin
    
    if isempty(idx1)
        error(['The energy value %g in the output group structure is not '...
                       'in the input group structure. The elements of the output'...
                       ' group structure must be a subset of the elements in the'...
                       ' input group structure.'],groupStructureOut(i))
    end
    
    if isempty(idx2)
        error(['The energy value %g in the output group structure is not'...
                       'in the input group structure. The elements of the output'...
                       ' group structure must be a subset of the elements in the'...
                       ' input group structure.'],groupStructureOut(i+1))
    end
    
    if ~isscalar(idx1)
        error(['The energy value %g in the output group structure is found '...
                        'multiple times in the input group structure'],groupStructureOut(i))
    end
    
    if ~isscalar(idx2)
        error(['The energy value %g in the output group structure is found'...
                        'multiple times in the input group structure'],groupStructureOut(i+1))
    end
    
    idx2 = idx2 - 1; % idx2 now gives the correct right edge of the i-th output bin
    
    % a more efficient way to do this would be to not walk through the
    % entire group structure each time; just look for the next value etc
    
    % assemble; give to GSiidxer{i} all the indices from idx1 to idx2
    GSIidxer{i} = idx1:idx2;
    
end

%% Collapse flux

% deal with flux first to save on computation
fluxOut = zeros(numBinsOut,1);

for i = 1:numBinsOut

    fluxOut(i) = sum( fluxIn( GSIidxer{i} ) );
    
end

if ~isstruct(xsIn) % execute normally % MUST MODIFY THIS FOR ERROR checking
    
    %% more error-checking
    
    % get the size of the input xs; either column vector or (square) kernel
    % matrix.
    [numRows_xsIn numCols_xsIn] = size(xsIn);
    
    if numRows_xsIn ~= numBinsIn % the xs must be the same length (size) as input bins
        error('The input group structure and cross section data do not have a similar size.')
    end
    
    if (numCols_xsIn ~=1) && (numCols_xsIn  ~= numRows_xsIn)
        % you're either a column vector or a kernel
        error('The input cross section data is not a kernel and is not a column vector.')
    end
    
    % error check on fluxIn
    [numRows_fluxIn numCols_fluxIn] = size(fluxIn);
    if (numRows_xsIn ~= numRows_fluxIn) || (numCols_fluxIn ~= 1)
        error(['The input flux does not have the same length as the input cross'...
            ' section data, or the input flux is not a column vector.'])
    end
    
    
    
    %% COLLAPSEIT!..the cross section data
    
    if numCols_xsIn == 1
        % vector collapsing
        
        % initialize
        xsOut = zeros(numBinsOut,1);
        
        for i = 1:numBinsOut % there's one xs per bin, not per edge
            idxs = GSIidxer{i};
            
            xsOut(i) = sum( xsIn(idxs) .* fluxIn(idxs) ) / fluxOut(i);
            
        end
        
    else
        % kernel collapsing
        
        % initialize
        xsOut = zeros(numBinsOut);
        
        for i = 1:numBinsOut % go through rows, E
            for j = i:numBinsOut % go through columns, E'
                
                rowidxs = GSIidxer{i}; % indexes rows
                colidxs = GSIidxer{j}; % indexes columns; we pull this flux
                
                % extract portion of xsIn kernel, row-collapse, and make into column vector
                xsIn_rowcollapsed = sum( xsIn(rowidxs,colidxs) )';
                
                xsOut(i,j) = sum( xsIn_rowcollapsed .* fluxIn(colidxs)) / fluxOut(j);
                
            end
        end
    end

% STRUCTURE
elseif isstruct(xsIn)

%     xstypes = {'scatker','fission','absorb','transport','vfission'};
    % xsOut is a struct that we create with these 4 fields, size out will be size in
    xsOut = struct(xstypes{1},[],xstypes{2},[],xstypes{3},[],xstypes{4},[]);
    
    numXSTypes = length(xstypes); % we do this service for each of these four xs's

    % fluxIn is a column vector
    
    for k = 1:numXSTypes
        
        thisXS = xsIn.(xstypes{k});

        [numRows_thisXS, numCols_thisXS] = size(thisXS);

        if numCols_thisXS == 1
            % vector collapsing

            xsOut.(xstypes{k}) = zeros(numBinsOut,1);
            
            for i = 1:numBinsOut
                idxs = GSIidxer{i};
            
                xsOut.(xstypes{k})(i) = sum( thisXS(idxs) .* fluxIn(idxs) ) / fluxOut(i);
            
            end
            
        else
            % kernel collapsing
            
            xsOut.(xstypes{k}) = zeros(numBinsOut);
            
            for i = 1:numBinsOut % go through rows, E
                
                for j = i:numBinsOut % go through columns, E'
                
                    rowidxs = GSIidxer{i}; % indexes rows
                    colidxs = GSIidxer{j}; % indexes columns; we pull this flux
                    
                    % extract portion of xsIn kernel, row-collapse, andmake into column vector
                    xsIn_rowcollapsed = sum( thisXS(rowidxs,colidxs) )';
                    
                    xsOut.(xstypes{k})(i,j) = sum( xsIn_rowcollapsed .* fluxIn(colidxs)) / fluxOut(j);
                
                end
            end
        end
    end

else
    error('The input cross section data is not in the correct format.') % improve this guy
end


end % end of function


%% TRASH
