function Results = Diffusion2D(DiffParamIn,g2D,g,p)
%DIFFUSION2D Creates mapping of regions to grid points for 2D
%diffusion, for an assembly-level reactor geometry.
% Description:  From the geometry structure g and parameter structure p,
% and the internal definition of assemblies and regions to assemblies,
% this function provides the physical dimensions of the reactor core,
% defines the numerical grid, and assigns region-defined material
% properties to grid points. The idea is that createMap4 is run to generate
% a specific reactor geometry. Other functions would be used to generate
% slightly different geometries, such as ones with different radii or four
% core regions. g2D (geometry 2D) structure has fields
%
% g2D.
%    .regionMap     holds the region assignment of each grid point
%    .edgePhys      side length (cm) of reactor, including reflector
%    .assy          fairly detailed structure defining all assemblies
%
% USE:  Results = Diffusion2D(DiffParamIn,g2D,g,p)
%
% NOTES: It is intended to use an improved and corrected function
% createMap5, which ensures that there is a uniform number of grid points
% in each assembly. Side length of each assembly is hard-coded as 22 cm.
% There are 180 assemblies, with an equal number of assemblies in each of 3
% regions. Region-assembly assignments can be modified on an individual
% basis, i.e. a water assembly can be placed in the center of the core.
% Assemblies are butted up against each other.
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110520 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%   1- finish implementing richer boundary conditions
%
% DEPENDENCIES: none
%   drawAssy (subfunction)
%   fillAssy (subfunction)
%

disp('Now Diffusing...');
global numNodeI numNodeJ numIntNodeI numIntNodeJ

numGroups = p.nGroups;

regionMap = g2D.regionMap;
edgePhys = g2D.edgePhys;

% start timing
tic

%% GEOMETRY INPUT
disp('Now arranging geometry and parameters...');

% this variable controls what is plotted 
% (1) enables any plots
% (2) plot dots at the grid points
% (3) plot lines connecting grid points
% (4) label grid points with their p-coordinate (linear coordinate)
% (5) label regions/campaigns with the fuel type or moderator (1-4) at that
% grid point
% (6)fill square around grid points according to region/campaign there
plotCore = [0 0 0 0 0 0];

%% BOUNDARY CONDITIONS
BoundCond = [0 0 0 0]; % left bottom right top. (0) zero flux, (1) symmetry

if sum(BoundCond > 1) > 0
    error(['Specified boundary conditions are incorrect. '...
        'Must be either 0 or 1 for each boundary']);
end

% initialize loop parameters / geometry
[numNodeI numNodeJ] = size(regionMap); % number of grid points in each direction
numNodes = numNodeI*numNodeJ;

% I IS THE Y-DIRECTION
% J IS THE X-DIRECTION

numIntNodeI = numNodeI - 2; % internal grid points
numIntNodeJ = numNodeJ - 2;
numIntNodes = numIntNodeI*numIntNodeJ;

%% Geometry

edgeSize = edgePhys; % cm
reactorLength = edgeSize; % cm
reactorWidth = edgeSize;
DeltaX = reactorWidth/(numNodeJ-1); % grid spacing
DeltaY = reactorLength/(numNodeI-1);

X = linspace(0,reactorWidth,numNodeJ); % physical coordinates of grid pts
Y = linspace(0,reactorLength,numNodeI);

% Plot the domain
close all
for k = 1
if plotCore(1)
    hGeom = figure(1);
    gGeom = axes;
    axis equal xy;
    set(gGeom,'XLim',[-reactorWidth*0.10 reactorWidth*1.10],'YLim',[-reactorLength*.10 reactorLength*1.10]);
    % can use axis ij and axis xy commands here maybe
    hold on;
    
    if plotCore(2)
        % plot dots for grid points
        for i = 1:numNodeI
            for j = 1:numNodeJ
                plot(X(j),Y(i),'k.','MarkerSize',15)
            end
        end
    end
    
    if plotCore(3)
        % plot lines for grid spaces
        offsetX = DeltaX * 0.15;
        offsetY = DeltaY * 0.15;
        % horizontal lines
        for i = 1:numNodeI
            for j = 1:numNodeJ-1
                plot([X(j)+offsetX, X(j+1)-offsetX],Y(i)*[1 1],'k')
            end
        end
    

        % vertical lines
        for i = 1:numNodeI-1
            for j = 1:numNodeJ
                plot(X(j)*[1 1],[Y(i)+offsetY, Y(i+1)-offsetY],'k')
            end
        end
    end
    
    if plotCore(4)
        % label grid points
        p = 1;
        for j = 2:numNodeJ-1
            for i = 2:numNodeI-1
                text(X(j),Y(i),sprintf('p=%d ',p),...
                    'HorizontalAlignment','right','VerticalAlignment','top')
                p = p + 1;
            end
        end
    end
    
    if plotCore(5)
        % label regions/campaigns
        for j = 1:numNodeJ
            for i = 1:numNodeI
                text(X(j),Y(i),sprintf(' %d',regionMap(i,j)),'Color',[1 0 0],...
                    'HorizontalAlignment','left','VerticalAlignment','bottom')
            end
        end
    end
    
    if plotCore(6)
        % fill regions/campaigns
        cIncr = 1/(max(max(regionMap))-min(min(regionMap))+1);
        for j = 1:numNodeJ
            for i = 1:numNodeI
                fill(X(j)+DeltaX*[-.5 .5 .5 -.5],Y(i)+DeltaY*[-.5 -.5 .5 .5],regionMap(i,j)*cIncr*[1 1 1])
            end
        end
        
    end
end
end

%% Process DiffParamIn
%    DiffParam{regidx}.scatker
%                     .xssink
%                     .vfission
%                     .diffco
%                     .fission
%                     .difflength

DiffParam = DiffParamIn;

% initialize diff params at internal grid points. these properties are not needed
% at the edge of the domain. each grid point has an energy dimension as
% well
ScatkerInt = zeros(numIntNodeI,numIntNodeJ,numGroups,numGroups);
VfissionInt = zeros(numIntNodeI,numIntNodeJ,numGroups);
FissionInt = zeros(numIntNodeI,numIntNodeJ,numGroups);

% map properties from campaigns to grid points
for j = 1:numIntNodeJ
    for i = 1:numIntNodeI
        campidx = regionMap(i+1,j+1); % campaign index
        % each interior grid point has its own scattering kernel!
        ScatkerInt(i,j,:,:) = DiffParam{campidx}.scatker;
        VfissionInt(i,j,:) = DiffParam{campidx}.vfission;
        FissionInt(i,j,:) = DiffParam{campidx}.fission;
    end
end

% these diff params are important at the edges
Xssink = zeros(numNodeI,numNodeJ,numGroups);
Diffco = zeros(numNodeI,numNodeJ,numGroups);
DiffLength = zeros(numNodeI,numNodeJ,numGroups);

% same old same old
for j = 1:numNodeJ
    for i = 1:numNodeI
        campidx = regionMap(i,j);
        
        Xssink(i,j,:) = DiffParam{campidx}.xssink;
        Diffco(i,j,:) = DiffParam{campidx}.diffco;
        DiffLength(i,j,:) = DiffParam{campidx}.difflength;
    end
end

maxGridSpacing = max( edgePhys/p.gp );
minGridSpacing = min( edgePhys/p.gp );

worstGridCriterion = maxGridSpacing / min(min(min(DiffLength)));

% put properties in a useful form for diffusion (stacey form)
A = zeros(numIntNodeI,numIntNodeJ,5,numGroups);
for groupidx = 1:numGroups
    
    D = Diffco(:,:,groupidx); % this line might be memory inefficient; i do it to simplify the look of the calculation lines within the loop.
    
    for j = 2:numNodeJ - 1
        for i = 2:numNodeI - 1
            
            % must do some boundary condition stuff here
            I = i-1; % internal indices
            J = j-1;
            
            % A has 4 dimensions: 2 for the 2D geometry, one for the type
            % of coefficient, and one for the energy group.
            
            A(I,J,1,groupidx) = Xssink(i,j,groupidx)...
                + 1/DeltaX^2 * ( mean( D(i,[j-1 j+1]) ) + D(i,j) )...
                + 1/DeltaY^2 * ( mean( D([i-1 i+1],j) ) + D(i,j) );
            
            A(I,J,2,groupidx) = -1/DeltaX^2 * mean( D(i-1:i,j) );
            
            A(I,J,3,groupidx) = -1/DeltaX^2 * mean( D(i:i+1,j) );
            
            A(I,J,4,groupidx) = -1/DeltaY^2 * mean( D(i,j-1:j) );
            
            A(I,J,5,groupidx) = -1/DeltaY^2 * mean( D(i,j:j+1) );
            
        end
    end
end
% can probably move the 1/Delta^2 outside of the loop if iso-grid spacing

% process boundary conditions
% it makes the most sense to do the above loop from 3: ()-2 and then
% process zero boundary condition here as well. i mean if the boundcond is
% zero, then A should be zero there anywayyy?
% this currently does not work
if BoundCond(1) == 1 % left
    
    j = 2;
    A(:,j-1,4,:) = zeros(numIntNodeI,1,1,numGroups);
    
    for groupidx = 1:numGroups
        
        D = Diffco(:,j:j+1,groupidx);
        
        for i = 2:numNodeI-1
            A(i-1,j-1,1,groupidx) = Xssink(i,j,groupidx)...
                + 1/DeltaX^2 * sum( D(i,:) )...
                + 1/DeltaY^2 * ( mean( D([i-1 i+1],1) ) + D(i,1) );
        end
    end
end

if BoundCond(2) == 1 % bottom
    
    i = 2;
    A(i-1,:,2,:) = zeros(1,numIntNodeJ,1,numGroups);
    
    for groupidx = 1:numGroups
        
        D = Diffco(i:i+1,:,groupidx);
        
        for j = 2:numNodeJ-1
            A(i-1,j-1,1,groupidx) = Xssink(i,j,groupidx)...
                + 1/DeltaX^2 * ( mean( D(1,[j-1 j+1]) ) + D(1,j) )...
                + 1/DeltaY^2 * sum( D(:,j) );
        end
    end
end

if BoundCond(3) == 1 % right
    
    j = numNodeJ-1;
    A(:,j-1,5,:) = zeros(numIntNodeI,1,1,numGroups);
    
    for groupidx = 1:numGroups
        
        D = Diffco(:,j-1:j,groupidx);
        
        for i = 2:numNodeI-1
            A(i-1,j-1,1,groupidx) = Xssink(i,j,groupidx)...
                + 1/DeltaX^2 * sum( D(i,:) )...
                + 1/DeltaY^2 * ( mean( D([i-1 i+1],2) ) + D(i,2) );
        end
    end
end

if BoundCond(4) == 1 % top
    
    i = numNodeI-1;
    A(i-1,:,3,:) = zeros(1,numIntNodeJ,1,numGroups);
    
    for groupidx = 1:numGroups
        
        D = Diffco(i-1:i,:,groupidx);
        
        for j = 2:numNodeJ-1
            A(i-1,j-1,1,groupidx) = Xssink(i,j,groupidx)...
                + 1/DeltaX^2 * ( mean( D(2,[j-1 j+1]) ) + D(2,j) )...
                + 1/DeltaY^2 * sum( D(:,j) );
        end
    end
end
% small computational inefficiencies in indexing here.

%% Energy Iterate

% create fast flux guess and corresponding fast source guess, only at
% interior nodes.

% 2D domain
[Xmesh, Ymesh] = meshgrid(X(2:end-1),Y(2:end-1));

% calculate the distance of each grid point to the reactor center
R = sqrt( (reactorWidth/2-Xmesh).^2 +(reactorLength/2-Ymesh).^2 );

% initialize 3-D internal flux (2 directional, 1 energy) (flux at boundary
% is determined by boundary conditions)
FluxInt = zeros(numIntNodeI,numIntNodeJ,numGroups);

for groupidx = 1:numGroups % initialize each energy group's flux
    FluxInt(:,:,groupidx) = cos(R/max(max(R)));
end
groupidx = 1; % fast flux
sourceFastInt = VfissionInt(:,:,groupidx).*FluxInt(:,:,groupidx);
sourceFastInt = sourceFastInt/max(max(sourceFastInt)); % normalize. are these two numbers small that we get roundoff?
    
% Prepare iteration
disp('Now starting energy iterations...')
threshold = 1e-5; % on error
isErrorful = 1; % continue condition
maxIterations = 10000;
counter = 1; % iteration counter
ErrrVector = []; % used to plot error as a function of iteration
while counter <= maxIterations && isErrorful
    
    %% FAST diffusion!
    groupidx = 1;
    % solve diffusion for fast group neutrons
    FluxInt(:,:,groupidx) = SingleGroup2Dout( sourceFastInt, FluxInt(:,:,groupidx),...
       A(:,:,:,groupidx), BoundCond);

    if p.nGroups == 3
    %% RESONANCE diffusion!
    groupidx = 2; % (from-group)
    gidx = 1; % (to-group)
    sourceResInt = zeros(numIntNodeI,numIntNodeJ,1);
    
    sourceResInt = ScatkerInt(:,:,gidx,groupidx).*FluxInt(:,:,gidx);

    FluxInt(:,:,groupidx) = SingleGroup2Dout( sourceResInt, FluxInt(:,:,groupidx),...
        A(:,:,:,groupidx), BoundCond);

    %% THERMAL Diffusion!
    groupidx = 3;
    gidx1 = 1;
    gidx2 = 2;
    sourceTherInt = zeros(numIntNodeI,numIntNodeJ,1);
    
    % there is some confusing matlab matrix magic here with 3rd-dimension
    % matrix products. it should take some time to understand how this line
    % adds up across energy groups
    sourceTherInt = ScatkerInt(:,:,gidx1,groupidx).*FluxInt(:,:,gidx1)...
                  + ScatkerInt(:,:,gidx2,groupidx).*FluxInt(:,:,gidx2);
    
    FluxInt(:,:,groupidx) = SingleGroup2Dout( sourceTherInt, FluxInt(:,:,groupidx),...
        A(:,:,:,groupidx), BoundCond);

    end
    %% Compute fast sources
    
    sourceFastIntNew = zeros(numIntNodeI,numIntNodeJ,1);
    
    for groupidx = 1:numGroups
        sourceFastIntNew = sourceFastIntNew + VfissionInt(:,:,groupidx).*FluxInt(:,:,groupidx);
    end
    
    Keff = norm(sourceFastIntNew)/norm(sourceFastInt);
   
    % IS THIS LINE NOT VALID ANMORE?
    sourceFastIntNew = sourceFastIntNew/max(max(sourceFastIntNew));

    errr = norm(sourceFastIntNew-sourceFastInt)/norm(sourceFastIntNew);
    
    ErrrVector = [ErrrVector errr]; % for plotting error

    if mod(counter,50) == 0
        fprintf('Iteration %d\n',counter);
    end
    
    % clean up and start over
    counter = counter+1;
    isErrorful = errr > threshold;
    sourceFastInt = sourceFastIntNew; 
    
end

%% CONVERGENCE
disp('Converged! Now plotting results...');

runtime = toc;

% output flux is for the internal nodes
fluxOut = FluxInt;
        
% effective power
powerEff = zeros(numIntNodeI,numIntNodeJ,numGroups);

powerEff = FissionInt.*FluxInt; % ASSUMES CONSTANT GRID SPACING.

powerEff = sum(powerEff,3);

powerTotal = sum(sum(powerEff));

powerFrac = zeros(g.nRegions,1);

for i = 1:g.nRegions
    powerFrac(i) = sum( sum( powerEff .* (regionMap(2:end-1,2:end-1) == i) ) ) / powerTotal;
end

% create Results
Results = struct('Xmesh',Xmesh,...
                 'Ymesh',Ymesh,...
        'minGridSpacing',minGridSpacing,...
        'maxGridSpacing',maxGridSpacing,...
    'worstGridCriterion',worstGridCriterion,...
                  'Flux',fluxOut,...
             'powerFrac',powerFrac,...
                  'keff',Keff,...
             'powerGrid',powerEff,...
               'runtime',runtime,...
           'nIterations',counter-1,...
                 'Error',ErrrVector,...
               'thetime',datestr(now,31));

end % ends everything


% the function below is no longer used
% the function below is no longer used
% the function below is no longer used
function [fluxIntOut] = SingleGroup2D(source,fluxIntA,A,BoundCond) 

% THIS FUNCTION IS DEPRACATED. I NO LONGER USE THIS SUBFUNCTION;
% SINGLEGROUP2DOUT HAS REPLACED THIS

global numNodeI numNodeJ numIntNodeI numIntNodeJ % A

% global A
% Calculate A, the coefficient matrix
% A = zeros(numSpaces,numSpaces);

% if too many iterations, do something bad.

%% Process Boundary conditions

fluxA = zeros(numNodeI,numNodeJ);

% these assignments are wrong, and assume a symmetry condition.
% for bcidx = 1:4
% if BoundCond(bcidx) == 0
% else if. % error check to make sure only 0's or 1's are passed?

% zero
% left
if BoundCond(1) == 0
    fluxA(1:numNodeI,1) = zeros(numNodeI,1);
end

% bottom (physically)
if BoundCond(2) == 0
    fluxA(1,1:numNodeJ) = zeros(1,numNodeJ);
end

% right
if BoundCond(3) == 0
    fluxA(1:numNodeI,numNodeJ) = zeros(numNodeI,1);
end

% top (physically)
if BoundCond(4) == 0
    fluxA(numNodeI,1:numNodeJ) = zeros(1,numNodeJ);
end

% symmetry
if BoundCond(1) == 1
    fluxA(2:numNodeI-1,1) = fluxA(2:numNodeI-1,2);
    fluxA(1,1) = fluxA(2,2);
    fluxA(numNodeI,1) = fluxA(numNodeI-1,2);
end

if BoundCond(2) == 1
    fluxA(1,2:numNodeJ-1) = fluxA(2,2:numNodeJ-1);
    fluxA(1,1) = fluxA(2,2);
    fluxA(1,numNodeJ) = fluxA(2,numNodeJ-1);    
end

if BoundCond(3) == 1
    fluxA(2:numNodeI-1,numNodeJ) = fluxA(2:numNodeI-1,numNodeJ-1);
    fluxA(1,1) = fluxA(2,2);
    fluxA(numNodeI,numNodeJ) = fluxA(numNodeI-1,numNodeJ-1);
end

if BoundCond(4) == 1
    fluxA(numNodeI,2:numNodeJ-1) = fluxA(numNodeI-1,2:numNodeJ-1);
    fluxA(numNodeI,1) = fluxA(numNodeI-1,2);
    fluxA(numNodeI,numNodeJ) = fluxA(numNodeI-1,numNodeJ-1);
end
% THE SYMMETRY CONDITION MUST BE CONTINUALLY UPDATED INSIDE THE LOOP
fluxB = fluxA;

% set up iteration
maxIterations = 2000;
isErrorful = 1;
itcounter = 0;
thresholdp = 1e-2;
while isErrorful

    p = 1;
    for j = 2:numIntNodeJ
        for i = 2:numIntNodeI
            
            % i'm not so sure this is correct with a and b. is there any
            % overlap/overwiting? there would be overwriting if i used a p
            % notation. i should look up the gauss-seidel / successive
            % relaxation method.
            % this calculation requires knowledge of boundary fluxes.
            fluxB(i,j) = source(i,j)...
                       - A(p,2)*fluxA(i-1,j) - A(p,4)*fluxA(i,j-1)...
                       - A(p,3)*fluxA(i+1,j) - A(p,5)*fluxA(i,j+1);
                          
            fluxB(i,j) = fluxB(i,j)/A(p,1);
            
            p = p + 1;
            
            % i need to change A according to boundary conditions
            
        end
    end
    
    % take care of boundary conditions; assign to boundary the correct flux.
    
    isErrorful = sum(sum( abs(fluxA - fluxB)./fluxA > thresholdp )) > 0 ;
%     isErrorful = norm( (fluxA-fluxB)./fluxA ) > thresholdp;

    if itcounter > maxIterations
        error('SingleGroup2D has not converged. put more details in this error msg')
    end
    
    itcounter = itcounter + 1;
    fluxA = fluxB;
    
%     source = vfiss.*fluxIntA; % should this ever be updated inside? is the source considered known?
    
end

fluxIntOut = fluxB(2:numNodeI-1,2:numNodeJ-1);

end
