function g2D = createMapSimple(g,p)
%CREATEMAPSIMPLE Creates mapping of regions to grid points for 2D
%diffusion, for a simple concentric-circles reactor geometry.
% Description:  From the geometry structure g and parameter structure p,
% this function provides the physical dimensions of the reactor core,
% defines the numerical grid, and assigns region-defined material
% properties to grid points. g2D (geometry 2D) structure has fields
%
% g2D.
%    .regionMap     holds the region assignment of each grid point
%    .edgePhys      side length (cm) of reactor, including reflector
%
% USE:  g2D = createMapSimple(g,p)
%
% NOTES: Physical size of reflector and of reactor core are hard-coded.
% Side length of reactor is 700 cm and core radius is 200 cm. Grid spacing
% is uniform and isotropic.
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110520 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%
% DEPENDENCIES: none
%

% plot region map array as a surface plot
doPlot = 0;

edgePhys = 700; %cm
edgeNode = p.gp;

cmPerNode = edgePhys/edgeNode;

% initialize regionMap as reflector region
regionMap = 4*ones(round(edgeNode));

% radius of core is hard-coded
coreRadius = 200; % cm

% radiuf campaign (region) boundaries
campRadius = sqrt(coreRadius^2/3*(3:-1:1));

% physical coordinate (in both x and y) of center of core
center = edgePhys/2;

% assign region indices to regionMap by approximating circles on the
% rectangular grid.
for regidx = 1:3
    for y = -campRadius(regidx):campRadius(regidx)
        X = sqrt(campRadius(regidx)^2-y^2);
        for x = -X:X
            
            % this rudimentary calculation exploits the fact that 
            i = y+center;
            j = x+center;
            i = round(i/cmPerNode);
            j = round(j/cmPerNode);
            regionMap(i,j) = regidx;
        end
    end
end

close all

if doPlot
    figure;
    colormap(summer);

    surf(regionMap,'EdgeColor','none');
    axis equal;
    view(0,90);
end

% create output structure
g2D = struct;
g2D.regionMap = regionMap;
g2D.edgePhys = edgePhys;

end