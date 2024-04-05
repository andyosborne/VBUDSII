function f = plot2DSlices(h,g,g2,Results,p,nSlices,varargin)
%PLOT2DSLICES Takes the results of 2D diffusion and plots radial flux
%slices.
% Description:  In order to make easy a comparison between 1D and 2D flux
% solutions, this function takes handles for a figure, an axes, and an
% inset axes, a structure Results, a parameter structure p (for nGroups),
% and some varargins and plots radial flux from the center (middle) of the
% 2D square grid. The directions in which the slices were taken are
% described by inset axes in the upper right of the figure. The first
% varargin is the style with which to plot the fluxes, e.g. 'k-'. This
% works for both 1-group and 3-group fluxes. Grouptoplot is for 3-group
% plots, and one out of the 3 groups can be selected for the plot. doInset
% controls whether the inset plot is plotted.
%    
% USE:  f = plot2DSlices(h,g,g2,Results,p,nSlices)
%       f = plot2DSlices(h,g,g2,Results,p,nSlices,style)
%       f = plot2DSlices(h,g,g2,Results,p,nSlices,style,grouptoplot)
%       f = plot2DSlices(h,g,g2,Results,p,nSlices,style,grouptoplot,doInset)
%
% NOTES: 
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110522 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%
% DEPENDENCIES: none
%

if isempty(varargin)
    % specify style
    % plot total flux
    Fluxor = Results.Flux;
    if p.nGroups == 1
        style = 'k-';
    elseif p.nGroups == 3
        style = {'r','g','b'};
    end
    % plot the inset figure
    doInset = 1;
elseif length(varargin) == 1
    % use specified style, which is hopefully consistent with p.nGroups
    Fluxor = Results.Flux;
    style = varargin{1};
    doInset = 1;
elseif length(varargin) == 2
    % varargin chooses which layer of the flux to plot
    style = varargin{1};
    Fluxor = Results.Flux(:,:,varargin{2});
    doInset = 1;
elseif length(varargin) == 3
    % varargin chooses if the inset plot is plotted
    style = varargin{1};
    Fluxor = Results.Flux(:,:,varargin{2});
    doInset = varargin{3};
end

% plotting domain
Xmesh = Results.Xmesh;
Ymesh = Results.Ymesh;

% width and length of reactor in physical units
Widthby2 = 1/2*Xmesh(1,end);
Lengthby2 = 1/2*Ymesh(end,1);

% unrotated direction of flux slice
Ray = [floor(.90*Widthby2); 0];

figure(h);
axes(g);
hold on; box on;

% choose the Angles at which the slices slice the flux
% Angles = pi/2*(0:3);
% Angles = pi/4*(0:7);
Angles = pi*[0 1/8 1/4 1/2 3/4];  % <<<======== choose what slices to plot!!!
Angles = Angles(1:min(nSlices,length(Angles)));
for idx = Angles
    
    % rotate the ray
    Ray2 = [cos(idx) -sin(idx); sin(idx) cos(idx)]*Ray; % rotations
    
    X = linspace(0,Ray2(1),100);
    Y = linspace(0,Ray2(2),100);
    
    R = sqrt(X.^2 + Y.^2);      % construct abscissa
    
    for groupidx = 1:p.nGroups  % interpolate
        Fluxthis = interp2(Xmesh,Ymesh,Fluxor(:,:,groupidx),X+Widthby2,(Y+Lengthby2)');
        Fluxline = diag(Fluxthis);
        if p.nGroups == 3
            f = plot(R,Fluxline,style{groupidx}); %colorstr{groupidx});
        else
            f = plot(R,Fluxline,style);
        end

    end
end

% label the plot
ylabel('Flux \phi (cm^{-2} s^{-1})')
xlabel('R (cm)')
if p.nGroups == 3
    legend('fast','resonation','thermal','Location','NorthEast');
end

% plot inset figure showing what slices were plotted
axes(g2);
set(g2,'Position',[.70 .55 .15 .15],'Layer','top','Box','off',...
    'XTick',[],'YTick',[]);
axis equal off;
hold on;
if doInset
    plot(0,0,'ro','MarkerSize',10);
    for idx = Angles
        Ray = [cos(idx) -sin(idx); sin(idx) cos(idx)]*[0 1]';
        
        plot([0 Ray(1)],[0 Ray(2)],'k');
    end
end