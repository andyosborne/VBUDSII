function [Xmesh,Ymesh,Fluxor] = MapCreator5(DiffParam,DiffParamRefl,gp)

% last edited (by and date): cld72@cornell.edu M110213

global AssySide

numGroups = 3;
% function [regionMap,edgeSize] = MapCreator2

% THIS COMMENT BELOW: future of MapCreator UI, etc.
% creates a figure window and allows ginput of what cell goes where
% it knows how many cells/regions there are, and allows you to click where
% cell 1 is, then where cell 2 is, then where cell 3 is. then the
% reflector. it knows how fine we want the grid to be, and accordingly puts
% the right cell identifier to the correct clicked area of the reactor
% core.

doPlot = 0; % plot input 
doPlot2 = 0; % plot output

%%%%%%%%%%%%%%%
%% RXR SETUP %%
%%%%%%%%%%%%%%%

% plot reactor boundaries
gpPerAssySide = 4;
gpPerAssy = gpPerAssySide^2;
assyPerSide = 24;
gp = assyPerSide*gpPerAssySide;
AssySide = 22; % side length of assy (cm)

fprintf('Number of gridpoints per assembly side: %d\n',gpPerAssySide);
fprintf('Number of gridpoints per assembly: %d\n',gpPerAssy);
fprintf('Number of gridpoints along reactor side: %d\n',gp);
fprintf('Total number of gridpoints: %d\n',gp^2);

edgeNode = gp; % THIS IS WHERE GRIDPOINTS ARE DEFINED: # of gridpoints on a side

edgePhys = assyPerSide*AssySide; % physical, (cm) % physical simulation size is an integer multiple of assy size
ReactorWidth = edgePhys;
ReactorLength = edgePhys;

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
    % set(g1,'XLim',ReactorWidth*[-sidebuff
    % sidebuff+1],'YLim',ReactorLength*[-sidebuff sidebuff+1])
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
    % set(g2,'XLim',ReactorWidth*[-sidebuff sidebuff+1],'YLim',ReactorLength*[-sidebuff sidebuff+1])
    
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

%%%%%%%%%%%%%%%%%%%
%% RUN DIFFUSION %%
%%%%%%%%%%%%%%%%%%%

% the exciting (and longgg) part!!!
[Xmesh,Ymesh,Fluxor,PowerEff] = Diffuse2D( DiffParam, DiffParamRefl, regionMapIn, edgePhys);

%%%%%%%%%%%%%%%%%
%% POSTPROCESS %%
%%%%%%%%%%%%%%%%%

% calculate normalized average power from each assy (from all 3 groups)
% go through each assembly.
maxPower = -inf;
for p = 1:length(Assy)
    PowerGrid = sum(PowerEff(Assy(p).iIdxr,Assy(p).jIdxr,:),3);
    Power = mean(mean(PowerGrid)); % averaging!
    Assy(p).power = Power;
    
    if Power > maxPower
        maxPower = Power;
    end
end

if doPlot2
   
    % top-down plot of assy's, color indicates relative power level
    h4 = figure;
    g4 = axes;
    axis equal off;
    hold on;

    plot(ReactorWidth*[0 1 1 0 0],ReactorLength*[0 0 1 1 0],'k')
    
    for p = 1:length(Assy)
        fillAssy(h4,Assy(p),Assy(p).power/maxPower);
    end
    
    % 3D plot of assy power levels    
    h5 = figure;
    g5 = axes;
    axis equal off;
    hold on;
    view(3);
    
    fill(ReactorWidth*[0 1 1 0 0],ReactorLength*[0 0 1 1 0],[.8 .8 .8])
    
    for p = 1:length(Assy)
        fill3Assy(h5,Assy(p),Assy(p).power/maxPower);
    end
    
end
end

function drawAssy(fhandle,Assy)

global AssySide

figure(fhandle)
plot(Assy.X+AssySide*[0 1 1 0 0],Assy.Y+AssySide*[0 0 1 1 0],'k')
    
end

function fillAssy(fhandle,Assy,color)

global AssySide

figure(fhandle)
fill(Assy.X+AssySide*[0 1 1 0],Assy.Y+AssySide*[0 0 1 1],color*[1 1 1])
    
end

function fill3Assy(fhandle,Assy,powah)

global AssySide maxPowah

m = .2;
figure(fhandle)
% X = Assy.X + AssySide*[ [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 0 0 0]', [0 1 1 0 0]' ];
% Y = Assy.Y + AssySide*[ [0 0 0 0 0]', [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 1 1 0]' ];
% Z = 10*AssySide*powah*[ [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [1 1 1 1 1]' ];
% fill3(X,Y,Z,Assy.cellid/4*[1 1 1])


X = Assy.X+m*AssySide/2 + (1-m)*AssySide*[ [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 0 0 0]', [0 1 1 0 0]' ];
Y = Assy.Y+m*AssySide/2 + (1-m)*AssySide*[ [0 0 0 0 0]', [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 1 1 0]' ];
Z = [ [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [1 1 1 1 1]' ];
fill3(X,Y,10*AssySide*powah*Z,Assy.cellid/4*[1 1 1])

m = 0;
X = Assy.X+m*AssySide/2 + (1-m)*AssySide*[ [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 0 0 0]', [0 1 1 0 0]' ];
Y = Assy.Y+m*AssySide/2 + (1-m)*AssySide*[ [0 0 0 0 0]', [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 1 1 0]' ];
% plot3(X,Y,10*AssySide*Z,'Color',[.9 .9 .9]);

end