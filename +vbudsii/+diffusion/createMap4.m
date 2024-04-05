function g2D = createMap4(g,p)
%CREATEMAP4 Creates mapping of regions to grid points for 2D
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
% USE:  g2D = createMap4(g,p)
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
%   1- must fix assignment of regions to assemblies (indices were switched
%   by accident
%   2- uniform # of grid points per assembly
%   3- decide how to manage boundary conditions: phantom grid points on
%   outside?
%
% DEPENDENCIES:
%   drawAssy (subfunction)
%   fillAssy (subfunction)
%

global AssySide

numGroups = p.nGroups;

doPlot = 0; % plot input 
doPlot2 = 0; % plot output

%% RXR SETUP

% plot reactor boundaries
edgePhys = 500; % physical, (cm)
ReactorWidth = edgePhys;
ReactorLength = edgePhys;

edgeNode = p.gp; % THIS IS WHERE GRIDPOINTS ARE DEFINED: # of gridpoints on a side

AssySide = 22; % side length of assy (cm)

sidebuff = .1; % not used currently; for plotting

% plot frame
if doPlot
    close all
    
    h1 = figure(1);
    set(h1,'Position',[100 100 1000 500])
    g1 = subplot(1,2,1);
    set(g1,'Position',[.05 .05 .44 .95])
    axis equal off;
    hold on;
    
    plot(ReactorWidth*[0 1 1 0 0],ReactorLength*[0 0 1 1 0],'k')

end

%% create Assy

% assy holds information about all assemblies; their location (physical and numerical)
% campaign, and power level (calculated by Diffuse2D
Assy = struct('X',[],'Y',[],'jIdxr',[],'iIdxr',[],'cellid',[],'power',[]);

p = 1; % struct index
R = 7; % parameter for creating proper number of assy's
for j = 0:R
    x = edgePhys/2-AssySide-AssySide*j;
    
    % define location of assy's and assign CAMPAIGNS
    for i = 0:sqrt(R^2-j^2)
        y = edgePhys/2-AssySide-AssySide*i;
        Assy(p).X = x;
        Assy(p).Y = y;
        if sqrt(i^2 + j^2) <= 3.7
            Assy(p).cellid = 3; % <----------campaign 3
        elseif sqrt(i^2 + j^2) <= 5.5
            Assy(p).cellid = 2; % <----------campaign 2
        else
            Assy(p).cellid = 1; % <----------campaign 1
        end
        
        if doPlot % plotting option
            drawAssy(h1,Assy(p));
            text(x+AssySide/2,y+AssySide/2,sprintf('%d',p),...
                'VerticalAlignment','middle','HorizontalAlignment','center');
        end
        p = p + 1;
    end
    
end

% copy above into quadrant IV
numQtrAssy = length(Assy);
for q = 1:numQtrAssy
    Assy(p) = Assy(q);
    Assy(p).X = Assy(q).X + 2*(edgePhys/2-Assy(q).X) - AssySide;
    if doPlot
    drawAssy(h1,Assy(p));
    text(Assy(p).X+AssySide/2,Assy(p).Y+AssySide/2,sprintf('%d',p),...
        'VerticalAlignment','middle','HorizontalAlignment','center');
    end
    p = p + 1;
end

% copy to quadrant II
for q = 1:numQtrAssy
    Assy(p) = Assy(q);
    Assy(p).Y = Assy(q).Y + 2*(edgePhys/2-Assy(q).Y) - AssySide;
    if doPlot
    drawAssy(h1,Assy(p));
    text(Assy(p).X+AssySide/2,Assy(p).Y+AssySide/2,sprintf('%d',p),...
        'VerticalAlignment','middle','HorizontalAlignment','center');
    end
    p = p + 1;
end

% copy to quadrant I
for q = 1:numQtrAssy
    Assy(p) = Assy(q);
    Assy(p).X = Assy(q).X + 2*(edgePhys/2-Assy(q).X) - AssySide;
    Assy(p).Y = Assy(q).Y + 2*(edgePhys/2-Assy(q).Y) - AssySide;
    if doPlot
    drawAssy(h1,Assy(p)); % subfunction below
    text(Assy(p).X+AssySide/2,Assy(p).Y+AssySide/2,sprintf('%d',p),...
        'VerticalAlignment','middle','HorizontalAlignment','center');
    end
    p = p + 1;
end

%% MAKE ASSY EDITS HERE
% Assy(88).cellid=4; % can make changes to assy campaign assignments
% use preceding plots to know which assy to edit

% plot assy cell assignments
if doPlot
    g2 = subplot(1,2,2);
    set(g2,'Position',[.5 .05 .44 .95])
    axis equal off;
    hold on;
    
    plot(ReactorWidth*[0 1 1 0 0],ReactorLength*[0 0 1 1 0],'k')
    
    for p = 1:length(Assy)
        fillAssy(h1,Assy(p),Assy(p).cellid/4); % subfunction below
        text(Assy(p).X+AssySide/2,Assy(p).Y+AssySide/2,sprintf('%d',Assy(p).cellid),...
            'VerticalAlignment','middle','HorizontalAlignment','center')
    end

    % put txt labels here
end


%% regionMap
% this is the information that Diffuse2D needs about grid point parameters
regionMap = 4*ones(edgeNode);

cmPerNode = edgePhys/edgeNode; % phys-to-numerical conversion
fprintf('Distance between grid points: %.5f cm\n',cmPerNode);
% ^^^ important output for checking diffusion length grid spacing criterion

% grid spacing
SideNode = round(AssySide/cmPerNode);

for p = 1:length(Assy) % go through each assembly
    x = Assy(p).X;
    y = Assy(p).Y;
    
    % determine which grid points to put in which assemblies
    I = round(y/cmPerNode);
    J = round(x/cmPerNode);
    
    Assy(p).iIdxr = I:I+SideNode;
    Assy(p).jIdxr = J:J+SideNode;
    
    for i = I:I+SideNode
        for j = J:J+SideNode
            regionMap(i,j) = Assy(p).cellid;
        end
    end
end

% expand regionMap: the above takes care of interior grid points, the below
% copies data near the boundary to the boundary. this data is used only
% for employing boundary conditions and is not physical
regionMapIn = zeros(edgeNode+2);
regionMapIn(2:end-1,2:end-1) = regionMap;
% sides
regionMapIn(2:end-1,1) = regionMap(:,1);
regionMapIn(1,2:end-1) = regionMap(1,:);
regionMapIn(2:end-1,end) = regionMap(:,end);
regionMapIn(end,2:end-1) = regionMap(end,:);
% corners
regionMapIn(1,1) = regionMap(1,1);
regionMapIn(end,1) = regionMap(end,1);
regionMapIn(end,end) = regionMap(end,end);
regionMapIn(1,end) = regionMap(1,end);

if doPlot
    h3 = figure;
    colormap(gray)
    surf(regionMapIn)
    axis equal;
    view(0,90);
end

% create output
g2D = struct;
g2D.regionMap = regionMap;
g2D.edgePhys = edgePhys;
g2D.assy = Assy;

end

function drawAssy(fhandle,Assy)
% draws an assembly in 2D

global AssySide

figure(fhandle)
plot(Assy.X+AssySide*[0 1 1 0 0],Assy.Y+AssySide*[0 0 1 1 0],'k')
    
end

function fillAssy(fhandle,Assy,color)
% fills an assembly in 2D with a certain color
global AssySide

figure(fhandle)
fill(Assy.X+AssySide*[0 1 1 0],Assy.Y+AssySide*[0 0 1 1],color*[1 1 1])
    
end

function fill3Assy(fhandle,Assy,powah)

global AssySide 

m = .2; % buffer space for clear plotting
figure(fhandle)

X = Assy.X+m*AssySide/2 + (1-m)*AssySide*[ [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 0 0 0]', [0 1 1 0 0]' ];
Y = Assy.Y+m*AssySide/2 + (1-m)*AssySide*[ [0 0 0 0 0]', [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 1 1 0]' ];
Z = [ [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [1 1 1 1 1]' ];
fill3(X,Y,10*AssySide*powah*Z,Assy.cellid/4*[1 1 1])

m = 0;
X = Assy.X+m*AssySide/2 + (1-m)*AssySide*[ [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 0 0 0]', [0 1 1 0 0]' ];
Y = Assy.Y+m*AssySide/2 + (1-m)*AssySide*[ [0 0 0 0 0]', [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 1 1 0]' ];

end