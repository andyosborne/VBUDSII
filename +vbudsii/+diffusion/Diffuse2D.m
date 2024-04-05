function [Xmesh,Ymesh,fluxOut,powerEff] = Diffuse2D(DiffParamIn,DiffParamRefl,regionMap,edgePhys)
%
% INPUTS
% DiffParamIn      struct containing 3-group diffusion parameters for each
%                  campaign
% DiffParamRefl    struct containing 3-group diffusion parameters for water
%                  reflector
% regionMap        matrix (N x N) containing the campaign designation of
%                  each grid point. also defines grid points and mesh size
% edgePhys         physical dimension of reactor's edge (e.g. 200 cm)

% OUTPUTS
% Xmesh,Ymesh      output of meshgrid, holds physical coordinates in 
%                  matrices for plotting
% fluxOut          matrix (N-2) x (N-2) flux at the interior grid points, 3-group
% powerEff         effective power at each grid point, using fission cross
%                  section

% last edited (by and date): cld72@cornell.edu M110213


disp('Now Diffusing...');
global numNodeI numNodeJ numIntNodeI numIntNodeJ

numGroups = 3;

doPlot = 0;


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
% DiffParam{campidx}.
%                  .scatker
%                  .xssink
%                  .vfission
%                  .diffco
%                  .fission

DiffParam = DiffParamIn;
DiffParam{4} = DiffParamRefl;

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

% same old same old
for j = 1:numNodeJ
    for i = 1:numNodeI
        campidx = regionMap(i,j);
        
        Xssink(i,j,:) = DiffParam{campidx}.xssink;
        Diffco(i,j,:) = DiffParam{campidx}.diffco;
    end
end

% put properties in a useful form for diffusion (stacey form)
A = zeros(numIntNodeI,numIntNodeJ,5,numGroups);
for groupidx = 1:numGroups
    
    D = Diffco(:,:,groupidx); % this line might be memory inefficient; i do it to simplify the look of the calculation lines within the loop.
    
    for j = 2:numNodeJ - 1
        for i = 2:numNodeI - 1
            
            % must do some boundary condition stuff here
            I = i-1; % internal indices
            J = j-1;
            
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

fluxOut = FluxInt;
        
if doPlot
    % Error plot
    hErrr = figure(2);
    gErrr = axes;
    semilogy(ErrrVector)
    xlabel('iteration (-)')
    ylabel('Error (-)')
    
    % surface plot for each energy group
    % additional plots are handled by main2D (for 2D flux slices) and
    % MapCreator4.
    hFlux = figure(3);
    colormap(copper);
    
    title(['Number of iterations: ' num2str(counter) '. Keff: ' num2str(Keff)])
    
    for groupidx = 1:numGroups
        subplot(1,3,groupidx)
        surf(Xmesh,Ymesh,fluxOut(:,:,groupidx),'EdgeColor','none');
        axis equal xy;
        view(0,90);
        xlabel(sprintf('Group %d',groupidx))
    end
end
fprintf('keff is: %.4f\n',Keff);

% effective power
powerEff = zeros(numIntNodeI,numIntNodeJ,numGroups);

powerEff = FissionInt.*FluxInt;




%% Output flux

% show runtime
fprintf('Number of iterations: %d\n',counter-1);

if runtime < 60
    fprintf('Runtime: %.2f seconds\n',runtime)
else
    fprintf('Runtime: %d min and %.2f seconds\n',...
        floor(runtime/60),mod(runtime,60))  
end

% fprintf(' ************** gp = %d\n',numNodeI-2)
% 
% fid = fopen('convergence_0220.txt','a+');
% fprintf(fid,'%d %.4d %d %.4f\n',numNodeI-2,runtime,counter-1,max(max(max(fluxOut))));
% fclose(fid);

end % ends everything


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

%% usefull

%         set(h,'NextPlot','replacechildren')

% This code performs the following operations:
%
      

% this initial flux guess does not belong here. where do we use our
% input flux?

%% trash

%     if mod(counter,10)==1
%         figure(h3);
%         plot(R*[1 1 1],fluxOut)
%         %         set(h3,'Visible','off')
%         xlabel(['iteration number ' num2str(counter) ', error is ' num2str(errr)])
%         legend('fast','res','thermal')
%         pause(.2)
%         
%     end
% 	% plot fast source
%     figure(h3);
%     plot(mean([R(2:end) R(1:end-1)],2),sourceFast)
%     xlabel('radial distance to center of space (cm)')
%     ylabel('Fast source')
%     set(h3,'Visible','off')

% %% Calculate power fractions
% 
% powerSpace = zeros(N,1);
% 
% for spaceidx = 1:N % can 
%     
%     nodalidx = spaceidx + 1;
%     
%     for groupidx = 1:numGroups % can do away with the energy loop using matrix algebra
%         
%         powerSpace(spaceidx) = powerSpace(spaceidx) + Vfission(spaceidx,groupidx)...
%                               * mean(fluxOut(nodalidx-1:nodalidx,groupidx));
%         
%     end
%     
% end
% 
% powerTotal = sum(powerSpace);
% 
% powerFrac = zeros(numCampaigns,1);
% 
% campaign3idxs = find(campaignMap==3);
% campaign2idxs = find(campaignMap==2);
% campaign1idxs = find(campaignMap==1);
% powerFrac(3) = sum( powerSpace( campaign3idxs ) ) / powerTotal; % campaign III
% powerFrac(2) = sum( powerSpace( campaign2idxs ) ) / powerTotal; % campaign II
% powerFrac(1) = sum( powerSpace( campaign1idxs ) ) / powerTotal; % campaign I
% 
% disp('Power Fractions!')
% for i = numCampaigns:-1:1
%     fprintf('Campaign %d: %.5f \n',i,powerFrac(i))
% end

% Scatker = zeros(numSpaceI,numSpaceJ,numGroups,numGroups);
% Xssink = zeros(numSpaceI,numSpaceJ,numGroups);
% Vfission = zeros(numSpaceI,numSpaceJ,numGroups);
% Diffco = zeros(numSpaceI,numSpaceJ,numGroups);
%     
% % assign properties to each grid point based on the grid point's
% % region/campaign.
% for spaceidxI = 1:numSpaceI
%     for spaceidxJ = 1:numSpaceJ
%         campidx = regionMap(spaceidxI,spaceidxJ);
%            
%         % each grid point has its own scattering kernel!
%         Scatker(spaceidxI,spaceidxJ,:,:) = DiffParam{campidx}.scatker;
%         
%         XSsink(spaceidxI,spaceidxJ,:) = DiffParam{campidx}.xssink;
%         Vfission(spaceidxI,spaceidxJ,:) = DiffParam{campidx}.vfission;
%         Diffco(spaceidxI,spaceidxJ,:) = DiffParam{campidx}.diffco;
%     end
% end

% A(p,p) = xssink(i,j)...
%         + 1/DeltaX^2 * ( mean( D(i,[j-1 j+1]) ) + D(i,j) )...
%         + 1/DeltaY^2 * ( mean( D([i-1 i+1],j) ) + D(i,j) );
%     
% A(p,p-1) = -1/DeltaX^2 * mean( D(i,j-1:j) );
% 
% A(p,p+1) = -1/DeltaX^2 * mean( D(i,j:j+1) );
% 
% A(p,p-I) = -1/DeltaY^2 * mean( D(i-1:i,j) );
% 
% A(p,p+I) = -1/DeltaY^2 * mean( D(i:i+1,j) );
% 
%             j = floor(p/I);
%     i = ceil(p/j);
% 
%     isErrorful = 0;
%     p = 1;
%     while p <= numIntNodes && ~isErrorful
%         isErrorful = abs(fluxA - fluxB)/fluxA > thresholdp;
%         p = p + 1;
%     end

% for p = 1:numNodesInt
%     [r, c] = ind2sub([numNodeI-2 numNodeJ-2],p);
%     FluxOutInt(r,c) = FluxInt(p);
% end

% fluxA(2:numNodeI-1,2:numNodeJ-1) = fluxIntA;
% fluxA(2:numNodeI-1,1) = fluxA(2:numNodeI-1,2);
% fluxA(1,2:numNodeJ-1) = fluxA(2,2:numNodeJ-1);
% fluxA(2:numNodeI-1,numNodeJ) = fluxA(2:numNodeI-1,numNodeJ-1);
% fluxA(numNodeI,2:numNodeJ-1) = fluxA(numNodeI-1,2:numNodeJ-1);
% % deal with 4 corners % what are the corner fluxes?
% fluxA(1,1) = fluxA(2,2);
% fluxA(1,numNodeJ) = fluxA(2,numNodeJ-1);
% fluxA(numNodeI,1) = fluxA(numNodeI-1,2);
% fluxA(numNodeI,numNodeJ) = fluxA(numNodeI-1,numNodeJ-1); % this is the tricky important one.
% 
% % boundary conditions
% % which has precedence at the corners? symmetry has precedence.
% % it would be good to not have to process these boundary conditions within
% % here; it's inefficient.
% % this is fairly inefficient right here i think
% bcIdxer = {1:numNodeI,1;
%            1,1:numNodeJ;
%            1:numNodeI,numNodeJ;
%            numNodeI,1:numNodeJ};% this might be the wrong side of the reactor.
% for bcidx = 1:4
%     if BoundCond(bcidx) == 0
%         fluxA(bcIdxer{bcidx,1},bcIdxr{bcidx,2}) = zeros(1);
%         
%     end
% end
% for j = 1:numIntNodeJ
%     for i = 1:numIntNodeI % calculate radius to each point
%         R(i,j) = sqrt( (reactorWidth/2-X(j+1))^2 + (reactorLength/2-Y(i+1))^2 );
%         % no real point in doing this vectorized
%     end
% end
% hFluxtot = figure(4);
% colormap(copper);
% surf(Xmesh,Ymesh,sum(fluxOut,3),'EdgeColor','none');
% axis xy;
% view(0,90);

% regionMapInput = [ 4 4 1 1 1;
%                    4 4 1 1 1;
%                    4 4 1 1 1;
%                    4 4 4 4 4;
%                    4 4 4 4 4];

% fid = fopen('Diffusion2DOutput.txt','a');
% % 3-GROUP ONLY
% fprintf(fid,['\r~~~~~~~ Run ' datestr(now,30) ' ~~~~~~~~~~~~~~\r']);
% for i = 1:N+1
%     fprintf(fid,'%7.5f   %7.5f   %7.5f\r',fluxOut(i,:));
% end
% 
% fclose(fid);
