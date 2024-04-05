function [xsOut, fluxOut] = CollapseRegion(relVolume, varargin)

global xstypes;

% RegionCollapse collapses cross section data for multiple cells to create
% a homogenized region using relative volume information and cell flux
% information. The resulting region cross sections are a result of a
% volume-flux average. This function accounts for the layout of fuel
% assemblies (e.g. distribution of IMF and UOX).
% RegionCollapse knows no material composition or energy structure.
%
% The function has been developed so that it must be called after
% CellCollapse, but can called before or after EnergyCollapse.
%
% INPUT
% relVolume           column vector (-) of the relative volume of each input cell
%                     in the current region, in order. sum(relVolume)=1 not
%                     required, the function normalizes this itself.
% varargin            variable number of cross section - flux (# cm^-2 s^-1)
%                     pairs, each defined for a given cell in the current
%                     region. Cross sections and fluxes can be of any group
%                     struct, but cross sections are macroscopic (cm^-1).
%                     Cross sections may be arrays, or structures with
%                     fields that contain either scalar or vector cross
%                     sections. In the latter case, the fields of the
%                     struct in XSTYPES are all operated on.
%                     Flux is always a vector with the same group structure
%                     as the cross section input.
% 
% OUTPUT
% xsOut               for array cross section input, this is an array with
%                     the same energy group structure as the input.
%                     for struct input, this is a struct with only the
%                     fields in XSTYPES, where the fields are arrays
%                     (including kernels) of cross sections with the same
%                     energy group structure as the input.
% fluxOut             volume-averaged flux for the region such that 
%         xsOut.*fluxOut = 1/sum(relVolume(relVolume(1)*xsIn1.*fluxIn1 + ...
%                                          relVolume(2)*xsIn2.*fluxIn2 + ...)
%
% Last updated (by and date): cld72@cornell.edu W101020

%% Process inputs

% error check

if nargin < 3
    error('This function requires at least three inputs')
end % check for an even number for (nargin-1)?

% check to see that relVolume is 1 dimensional

% the length of relVolume is the number of cells being given to RegionCollapse
[numCells numCols_V] = size(relVolume);

if ( numCols_V ~= 1 )
    error('The relative volume vector must be a column vector.')
end

if ( numCells ~= (nargin - 1)/2 )
    error('The length of the input relative volume must be the number of input cross section/flux pairs')
end

% normalize volume
relVolume = relVolume/sum(relVolume);

% initialize
if ~isstruct(varargin{1})
    % array input
    
    % grab size of input cross section (might be a kernel)
    [numRows_var numCols_var] = size(varargin{1}); 
    
    Cellxs = cell(numCells);
    Cellflux = zeros(numRows_var,numCells);
    
    for i = 1:numCells
        
        % shuffle cell array inputs into local variables
        Cellxs{i} = varargin{2*i-1}; % cross section for each cell.
                                     % use cell array because cross section could
                                     % could be a kernel (scatker)
        Cellflux(:,i) = varargin{2*i}; % flux for each cell

        if ( i < numCells ) && sum(size(varargin{2*i-1}) ~= size(varargin{2*i+1})) % beit matrix, scalar, or vector
            error('All input cross sections must be of the same dimensions.')
        end
        
    end
    
    %% do work
    
    % we are flux averaging as well!
    
    numer = zeros(numRows_var, numCols_var); % temporary numerator
    fluxOut = zeros(numRows_var,1);          % initialize fluxOut
    xsOut = zeros(numRows_var, numCols_var); % initialize xsOut
    
    for i = 1:numCells
        % add each cell's contribution to the function's output
        
        % the use of this intermediate term is important for managing kernels
        xstimesflux = zeros(numRows_var,numCols_var);
        
        for j = 1:numCols_var
            
            % intermediate term, multiply each column of the input cell
            % cross section by the flux for that cell, by energy group
            xstimesflux(:,j) = Cellxs{i}(:,j) .* Cellflux(:,i);
            
        end
        
        % weight cross section-flux by volume
        numer = numer + relVolume(i) * xstimesflux;
        
        % weight output flux by volume
        fluxOut = fluxOut + relVolume(i) * Cellflux(:,i);
        
    end
    
    % divide by flux for each energy level
    
    for i = 1:numCols_var
        xsOut(:,i) = numer(:,i) ./ fluxOut;
    end
        
elseif isstruct(varargin{1}) % need to do a more formal check. now we expct cellarray
    
    % xstypes needs to be made global, or use container.Map
%     xstypes = {'scatker','fission','absorb','transport','vfission'};
    
    % initialize output cross section struct
    xsOut = struct(xstypes{1},[],xstypes{2},[],xstypes{3},[],xstypes{4},[]);
    
    numXSTypes = length(xstypes); % the number of cross section types to work on
    
    numRows_var = length(varargin{2}); % number of energy groups
    
    Cellflux = zeros(numRows_var,numCells);
    for i = 1:numCells
        
        % shuffle cell array input into local variable
        Cellflux(:,i) = varargin{2*i}; % array
        
        % error check: flux inputs are of all the same size (group structure)
        
        if ( i < numCells ) && (sum(size(varargin{2*i}) ~= size(varargin{2*i+2})))
            error('All input fluxes must be of the same dimensions.')
        end
    end
    
    % take care of flux now because it doesn't depend on xstype
    fluxOut = zeros(numRows_var,1);
    for i = 1:numCells
        
        fluxOut = fluxOut + relVolume(i) * Cellflux(:,i);
        
    end
    
    % we depend on XS type now
    for k = 1:numXSTypes
        
        [numRows_var numCols_var] = size( varargin{1}.(xstypes{k}) ); % assuming all have the same size
        thisTCellxs = cell(1,numCells);
        
        for i = 1:numCells
                
            % shuffle cell array input into local cell array, holding cross
            % section arrays or kernel
            thisTCellxs{i} = varargin{2*i-1}.(xstypes{k}); % struct
            
            % error check: cross section inputs are of all the same size (group structure / kernel)
            if ( i < numCells ) && (sum(size(varargin{2*i-1}.(xstypes{k})) ~= size(varargin{2*i+1}.(xstypes{k}))))
                error('All input cross sections for this xstype must be of the same dimensions.')
            end

        end

        % initialize internal variable and output
        numer = zeros(numRows_var, numCols_var); % temporary
        xsOut.(xstypes{k}) = zeros(numRows_var, numCols_var);
        
        for i = 1:numCells
            
            % intermediate term
            xstimesflux = zeros(numRows_var,numCols_var);
            
            for j = 1:numCols_var

                xstimesflux(:,j) = thisTCellxs{i}(:,j) .* Cellflux(:,i);
                
            end
            
            numer = numer + relVolume(i) * xstimesflux;
            
        end
        
        % divide by flux for each energy level
        for i = 1:numCols_var
            
            xsOut.(xstypes{k})(:,i) = numer(:,i) ./ fluxOut;
            
        end
        
    end

else
    error('The input cross section data is not in the correct format.')
end

end

%% TRASH
%         if size(Cellxs{i} ~= size(Cellflux(:,i)) && (size(varargin{2*i}) ~= size(varargin{2*(i+1)}))
%             error('All input cross sections and fluxes must be of the same dimensions.')
%         end
%         error('The input cross sections and fluxes must be scalars or column vectors.')
%     end
% 
%     if (numCells ~= length(varargin{i})) || (numCells ~= length(varargin{i+1}))
%         error('All input cross sections and fluxes must be data of the same dimensions.')
%     end
    

    
%     if ( numCols_var == 1 )
%         numer = zeros(numRows_var,1);
%         
%         for i = 1:numCells
%             %               scalar         vec or scalar  vec or scalar
%             numer = numer + relVolume(i) * Cellxs{i} .* Cellflux(:,i);
%             
%             fluxOut = fluxOut + relVolume(i) * Cellflux(:,i);
%             
%         end
%         
%         xsOut = numer/fluxOut;
    
%     elseif ( numCols_var == numRows_var )






        
        
        
        
        
%         if (size(varargin{2*i-1}) ~= size(varargin{2*i+1})) % beit matrix, scalar, or vector
%             error('All input cross sections must be of the same dimensions.')
%         end
%         
%             Cellxs = cell(numCells);
%     Cellflux = zeros(numRows_var,numCells);
%     
%     for i = 1:numCells-1 % -1 to not not over-index by 1
%         
%         % in here we would deal with a struct input
%         
%         % shuffle cell array inputs into local holders
%         Cellxs{i} = varargin{2*i-1};
%         Cellflux(:,i) = varargin{2*i};
%         
%         if (size(varargin{2*i-1}) ~= size(varargin{2*i+1})) % beit matrix, scalar, or vector
%             error('All input cross sections must be of the same dimensions.')
%         end
%         
%     end
%     
%     
%     
%     
%     
%     
%     
%     
%                 % vector collapsing!
%                 for i = 1:numCells - 1
%                     
%                     Cellxs{i} = varargin{2*i-1}.(xstypes{k});
%                     
%                     
%                 end
%                 
%                     for j = 1:numCols_var
%         
%         % intermediate term
%         intermediateTerm1(:,j) = Cellxs{i}(:,j) .* Cellflux(:,i);
%         
%     end
% 
%             
%     end
% 
% 
% 
% 
%     CellParts = varargin{1}; % U235, U238, O16, or H2O
%     
%     xstypes = {'scatker','fission','absorb','transport'};
%     % xsOut is a struct that we create with these 4 fields, size out will be size in
%     xsOut = struct(xstypes{1},[],xstypes{2},[],xstypes{3},[],xstypes{4},[]);
%     
%     numXSTypes = length(xstypes); % we do this service for each of these four xs's
%     
%     numParts = length(CellParts);
%     
%     for i = 1:numXSTypes
%         for j = 1:numParts
%             xOut.(xstypes{i}) = xOut.(xstypes{i}) + Composition(j) * CellParts{j}.(xstypes{i});
%             
%         end
%     end
%     


% if sum(relVolume) ~= 1
%     error('The sum of the relative volume vector is not 1.')
% end










