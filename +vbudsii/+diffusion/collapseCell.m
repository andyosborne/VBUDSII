function [xsOut] = collapseCell(Composition,varargin)

global xstypes;

% CellCollapse collapses cross section data for multiple isotopes (unit
% substances) to a single molecule (material) using Composition
% information. This function thus accounts for the changing enrichment of
% the fuel.
% CellCollapse knows no geometry or energy group structure.
% 
% The function has been developed so as to be flexible with application,
% but in its intended use, the inputs are 110-group microscopic cross
% sections organized in a cell array of structs.
%
% Purpose:
%
% INPUT
% Composition         vector (# b^-1 cm^-1) of atomic densities corresponding to
%                     the input microscopic cross section data. generated
%                     by Bateman
% varargin            this function takes variable arguments of either
%                     multiple arrays (one for each entry of Composition)
%                     or a cell array of structs, where each cell's
%                     struct provides a list of cross sections to be
%                     collapsed
% multiple arrays     a scalar, vector, or kernel (barns) for each entry of
%                     Composition, with any (but consistent) energy group structure
% cell array          each cell holds a struct that corresponds to an entry
%                     of Composition, where the struct contains multiple
%                     cross sections. The cross sections that are operated
%                     on are given in the cell array XSTYPES
%
% OUTPUT
% xsOut               for multiple arrays input, this is an array (cm^-1) of the
%                     same size as the input cross sections, as macroscopic
%                     cross sections of composition Composition
%                     for cell array input, this is a single struct of the
%                     same size as the input cross sections, with the cross
%                     sections given in XSTYPES converted to macroscopic
%                     (cm^-1)  cross sections
%
% HOW IT WORKS
%
% Last updated (by and date): cld72@cornell.edu M101018

%% Process inputs

include_variables

if iscell(varargin{1}) && length(varargin) ~= 1
    error('The input cannot be multiple cell arrays. Either vectors or a single cell holding structures.')
end

if ~iscell(varargin{1}) % as originally programmed, the first varargin is not a cell array, so it's probably 
    
% error check

% check to see that Composition is 1 dimensional

% assume argins come in xs-flux pairs
[numParts numCols_C] = size(Composition);

if ( numCols_C ~= 1 )
    error('The composition vector must be a column vector.')
end

if ( numParts ~= nargin - 1 )
    error('The length of the input composition must be the number of input cross sections')
end

% initialize

[numRows_var numCols_var] = size(varargin{1}); % numRows_var is the number of energy bins

if ( numCols_var ~=1 ) && ( numCols_var ~= numRows_var )
    error('The input cross sections must be scalars or column vectors, or square kernels.')
end

Partxs = cell(1,numParts);

for i = 1:numParts % -1 to not not over-index by 1

    Partxs{i} = varargin{i};
    
    if i ~= numParts && (size(varargin{i}) ~= size(varargin{i+1}))
        error('All input cross sections must be of the same dimensions.')
    end
    
end

%% do work

xsOut = zeros(numRows_var,numCols_var);

for i = 1:numParts
    
    xsOut = xsOut + Composition(i) * Partxs{i};

end

% STRUCT INPUT
elseif iscell(varargin{1})
    
    CellParts = varargin{1}; % U235, U238, O16, or H2O
    
%     xstypes = {'scatker','fission','absorb','transport','vfission','gamma'};
    % xsOut is a struct that we create with these 4 fields, size out will be size in
    xsOut = struct(xstypes{1},[],xstypes{2},[],xstypes{3},[],xstypes{4},[]);
    
    numXSTypes = length(xstypes); % we do this service for each of these four xs's
    
    numParts = length(CellParts);
    
    for i = 1:numXSTypes
        
        [numRows_var, numCols_var] = size(CellParts{1}.(xstypes{i})); % ERROR-CHECK: this assumes consistency
        xsOut.(xstypes{i}) = zeros(numRows_var, numCols_var); % manages scatkers
        
        for j = 1:numParts

            xsOut.(xstypes{i}) = xsOut.(xstypes{i}) + Composition(j) * CellParts{j}.(xstypes{i});
            
        end
    end
    
    
    
% bad input
else
    error('The input cross section data is not in the correct format.') % improve this guy
end

end

%% TRASH
% format of input structure: microxs{matl,temp,beta}.ngamma(E)
% isotopes/moleculares are called cell components

% possibly output the composition (atomic density) for the cell?
% want modular portable codes: use cellcollapse to simply make something
% macroscopic
% Composition accounts for enrichment

%     for i = 1:numParts
%         
% %         thisPart CellParts{i};
%         
%         for j = 1:numXSTypes
%             thisXSType = thisPart.(xstypes{j});
%             
%             % collapse!
%             for k = 1:numParts
%                 xsOut = xsOut + Composition(i) * thisXSType
%             
%             end
%         end

% uo2flux = FLUX(219:-2:1,1);
% h2oflux = FLUX(219:-2:1,2);
% flux = (h2oflux + uo2flux)/2;
% [U235microscat flux] = EnergyCollapse(edges,[1e-4 1 100e3 10e6],O16_0_2_600,flux)
% U235microscat;
%     xssize = length(thisCellParts{1}.fission); %this needs to be more
%     rigorous
% [xsOut, fluxOut] =