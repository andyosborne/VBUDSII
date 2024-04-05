function Results = Diffusion1D(DiffParamIn,g,p)
%DIFFUSION1D Solves 1-dimensional diffusion given geometry and parameter
%structures
% Description:  Solves either 1-group or 3-group neutron diffusion for an
% infinite cylinder core from 3-region (plus reflector) diffusion
% parameters, and reactor geometry specified in g, and run parameters
% specified in p. Produces a structure Results with the fields
%
% Results.R
%        .minGridSpacing
%        .maxGridSpacing
%        .worstGridCriterion
%        .Flux
%        .keff
%        .powerFrac
%        .runtime
%        .nIterations
%        .Error
%        .thetime
%    
% USE:  Results = Diffusion1D(DiffParamIn,g,p)
%
% NOTES: Method based off Dr. Cady 1994 formulation. Boundary conditions
% are hard-coded. Number of campaigns is hard-coded as 3.
%
%  \                           \             \       \      \
%   |         campaign III      |  camp II    |   I   | refl | 
%  /                           /             /       /      /
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110519 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%
% DEPENDENCIES:
%   SingleGroup (subfunction)
%   thomas
%

% global variables for SingleGroup
global N R Delta

% 1-group or 3-groups
numGroups = p.nGroups;

% start timing
tic

% option to plot the grid before diffusion starts. hard-coded as no.
plotCore = 0;
       
%% GEOMETRY INPUT
numCampaigns = g.nRegions -1;

myMapcase = 2;
gp = p.gp;

switch myMapcase
    case 1 % Dr. Cady's preferred grid
        campaignMap = [ 3 3 2 2 1 1 4 4 4]; 
    case 2 % Equal number of grid points per region
        campaignMap = [ 3*ones(1,gp) 2*ones(1,gp) 1*ones(1,gp) 4*ones(1,gp) ];
end

N = length(campaignMap); % number of grid points
k = 0:N; % grid indexer

%% Process DiffParamIn

DiffParam = DiffParamIn;

% form of DiffParam is:
%    DiffParam{regidx}.scatker
%                     .xssink
%                     .vfission
%                     .diffco
%                     .fission
%                     .difflength

Scatker = zeros(N,numGroups,numGroups);
Xssink = zeros(N,numGroups);
Vfission = zeros(N,numGroups);
Fission = zeros(N,numGroups);
Diffco = zeros(N,numGroups);
DiffLength = zeros(N,numGroups);
    
% map diffusion parameters from regions to grid points
for spaceidx = 1:N
    
    campidx = campaignMap(spaceidx);
    
    Scatker(spaceidx,:,:) = DiffParam{campidx}.scatker;
    
    Xssink(spaceidx,:) = DiffParam{campidx}.xssink;
    Vfission(spaceidx,:) = DiffParam{campidx}.vfission;
    Fission(spaceidx,:) = DiffParam{campidx}.fission;
    Diffco(spaceidx,:) = DiffParam{campidx}.diffco;
    DiffLength(spaceidx,:) = DiffParam{campidx}.difflength;
    
end

% actual geometry
CoreRadius = 200; % cm

CoreSpaces = N-sum(campaignMap==4); % side length of fuel elements (cm)
reflR = 1.3*CoreRadius;
campR = sqrt(CoreRadius^2 / 3 * (0:3)); % radius of campaign (region) boundaries

% define radius of grid points in each campaign
campaign3R = linspace(campR(1),campR(2),sum(campaignMap==3)+1);
campaign2R = linspace(campR(2),campR(3),sum(campaignMap==2)+1);
campaign1R = linspace(campR(3),campR(4),sum(campaignMap==1)+1);
reflR = linspace(campR(4),reflR,sum(campaignMap==4)+1);
R = [campaign3R campaign2R(2:end) campaign1R(2:end) reflR(2:end)];
R = R';

% plot the core
close all
if plotCore
    h = figure(1);
    for i = 1:N
        hold on;
        thet = 0:.01:2*pi;
        x = R(i)*cos(thet);
        y = R(i)*sin(thet);
        plot(x,y,'Color',[0 .5 .25])
        axis equal
    end
end

% thickenss of each annular space,
Delta = diff(R);

% grid spacing information for Results
minGridSpacing = min( Delta );
maxGridSpacing = max( Delta );

worstGridCriterion = max( Delta ./ min(DiffLength,[],2) );


%% BOUNDARY CONDITIONS
betaBC = [1 0];
extCurrent = [0 0];

%% ITERATION
% first fast source guess
avgdiffco = 0;
for i = 1:numCampaigns
    avgdiffco = avgdiffco + DiffParam{i}.diffco(1);
end
avgdiffco = avgdiffco / numCampaigns;
d = 2.13 * avgdiffco;
sourceFast = cos( R(2:end) / ( R (end) + d ) ); % space property

% error
threshold = 1e-5;
isErrorful = 1;
maxIterations = 2000;
counter = 1;
ErrrVector = []; % used to plot error as a function of iteration
while isErrorful && counter <= maxIterations
    
    fluxOut = zeros(N+1,numGroups);
    
    %% Fast group diffusion
    groupidx = 1;
    
    % fast
    fluxOut(:,groupidx) = SingleGroup( sourceFast, Xssink(:,groupidx), Diffco(:,groupidx), betaBC, extCurrent );
    
    if p.nGroups == 3
    %% Resonance diffusion!
    groupidx = 2;
    
    sourceRes = zeros(N,1);
    
    for spaceidx = 1:N
        nodalidx = spaceidx + 1;
        
        sourceRes(spaceidx) = sourceRes(spaceidx) + Scatker(spaceidx,groupidx-1,groupidx)...
                            * mean(fluxOut(nodalidx-1:nodalidx,1));
    end

    fluxOut(:,groupidx) = SingleGroup( sourceRes, Xssink(:,groupidx), Diffco(:,groupidx), betaBC, extCurrent);

    %% tHERMAL Diffusion!
    groupidx = 3;
    % 3group only!!!
    sourceTher = zeros(N,1);
        
    for spaceidx = 1:N
        nodalidx = spaceidx + 1;
        
        for gidx = 1:groupidx-1
            
            sourceTher(spaceidx) = sourceTher(spaceidx) + Scatker(spaceidx,gidx,groupidx)...
                                * mean(fluxOut(nodalidx-1:nodalidx,gidx));
            
        end
    end
    
    fluxOut(:,groupidx) = SingleGroup( sourceTher, Xssink(:,groupidx), Diffco(:,groupidx), betaBC, extCurrent);

    end
    %% Compute fast sources
    
    sourceFastNew = zeros(N,1);
    
    
    for spaceidx = 1:N
        for groupidx = 1:numGroups
            
            nodalidx = spaceidx + 1;
            
            sourceFastNew(spaceidx) = sourceFastNew(spaceidx) + Vfission(spaceidx,groupidx)...
                * mean(fluxOut(nodalidx-1:nodalidx,groupidx));
            
        end
    end
    
    Keff = norm(sourceFastNew)/norm(sourceFast);
  
    % normalize the source
    sourceFastNew = sourceFastNew/sourceFastNew(1);

    % relative error
    errr = norm(sourceFastNew-sourceFast)/norm(sourceFastNew);
    
    % store iteration error for results
    ErrrVector = [ErrrVector errr];
    %% Plot fluxes
    counter = counter+1;
    isErrorful = errr > threshold;
    sourceFast = sourceFastNew;
    
end

% we have convergence
runtime = toc;

%% Calculate power fractions

% initialize the power generated in each grid space across all energy
% groups
powerSpace = zeros(N,1);

for spaceidx = 1:N
    
    nodalidx = spaceidx + 1;
    
    for groupidx = 1:numGroups % can do away with the energy loop using matrix algebra
        
        powerSpace(spaceidx) = powerSpace(spaceidx) + Fission(spaceidx,groupidx)...
                              * mean(fluxOut(nodalidx-1:nodalidx,groupidx)) * Delta(spaceidx);
        
    end
    
end

powerTotal = sum(powerSpace);

powerFrac = zeros(numCampaigns,1);

campaign3idxs = find(campaignMap==3);
campaign2idxs = find(campaignMap==2);
campaign1idxs = find(campaignMap==1);
powerFrac(3) = sum( powerSpace( campaign3idxs ) ) / powerTotal; % campaign III
powerFrac(2) = sum( powerSpace( campaign2idxs ) ) / powerTotal; % campaign II
powerFrac(1) = sum( powerSpace( campaign1idxs ) ) / powerTotal; % campaign I

% prepare output
Results = struct('R',R,...
    'minGridSpacing',minGridSpacing,...
    'maxGridSpacing',maxGridSpacing,...
'worstGridCriterion',worstGridCriterion,...
              'Flux',fluxOut,...
              'keff',Keff,...
         'powerFrac',powerFrac,...
           'runtime',runtime,...
       'nIterations',counter-1,...
             'Error',ErrrVector,...
           'thetime',datestr(now,31));  
   
end

function [fluxOut] = SingleGroup(sourceIn,xssink,diffco,betaBC,extCurrent)

global N R Delta 

% from Dr. Cady 1994

%% GRIDSPACE LOOP!

% e, d, k defined for space
e = zeros(N,1);
d = zeros(N,1);
t = zeros(N,1);

% a, c, b, sigma defined for grid points
a = zeros(N,1);
c = zeros(N,1);

for spaceidx = 1:N % loop through spaces, get a, b, c, sigma
    
    % defined such that nodalidx = 1 when k = 0 (centerline) and
    % nodalidx = N+1 when k = N (boundary)
    nodalidx = spaceidx + 1;
    
    e(spaceidx) = 0.5 * xssink(spaceidx) * Delta(spaceidx);
    d(spaceidx) = 2.0 * diffco(spaceidx) / Delta(spaceidx);
    % for cylinders. t = 1/2 for slabs, r(nodalidx-1)^2 / (r(nodalidx)^2 +
    % r(nodalidx-1)^2) for spheres
    
    p = 1;
    
    t(spaceidx) = R(nodalidx-1)^p / ( R(nodalidx)^p + R(nodalidx-1)^p ); % R is defined nodally
    
    a(spaceidx) = (1 - t(spaceidx)) * ( e(spaceidx) - d(spaceidx) );
    c(spaceidx) = t(spaceidx)*( e(spaceidx) - d(spaceidx) );
    
end

% deal with b and sigma now
b = zeros(N+1,1);
sigma = zeros(N+1,1);

% centerline
nodalidx = 1; % k = 0
spaceidx = 1; % k = 1

b(nodalidx) = t(spaceidx) * e(spaceidx) + ( 1 - t(spaceidx) ) * d(spaceidx)...
             + 0.5 * (1 - betaBC(1)) / (1 + betaBC(1));

sigma(nodalidx) = 0.5 * sourceIn(spaceidx)*Delta(spaceidx) + 2 * extCurrent(1) / (1 + betaBC(1));

% edge of reactor (end of reflector...)
nodalidx = N+1;
spaceidx = N;

b(nodalidx) = (1 - t(spaceidx))*e(spaceidx) +  (t(spaceidx) * d(spaceidx))...
             + 0.5 * (1 - betaBC(2)) / (1 + betaBC(2));

sigma(nodalidx) = 0.5 * sourceIn(spaceidx)*Delta(spaceidx) + 2 * extCurrent(2) / (1 + betaBC(2));

for nodalidx = 2:N % nodalidx = 1 (k = 0) and nodalidx = N+1 (k = N) taken care of
    
    spaceidx = nodalidx-1;
    
    b(nodalidx) = (1 - t(spaceidx))*e(spaceidx) +  (t(spaceidx) * d(spaceidx))...
        + (t(spaceidx+1) * e(spaceidx+1)) + (1 - t(spaceidx+1))*d(spaceidx+1);
    
    sigma(nodalidx) = 0.5 * (sourceIn(spaceidx)*Delta(spaceidx) + sourceIn(spaceidx+1)*Delta(spaceidx+1));
    % again, source is defined on grid point nodes, and Delta is defined on
    % spaces
        
end

%% Solve the matrix!

    fluxOut = thomas(b, a, c, sigma);

end
