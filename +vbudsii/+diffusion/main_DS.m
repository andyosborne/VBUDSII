function DiffParam = main_DS(varargin)
fprintf('\n!!!New Run!!! at %s\n',datestr(now,31))
close all;

if length(varargin) == 2
    g = varargin{1};
    p = varargin{2};
elseif length(varargin) == 1 || isempty(varargin)
    if length(varargin) == 1
        p = varargin{1};
    elseif isempty(varargin)
        % DEFINE PARAMETERS
        p = struct('nDim',1,...
                'nGroups',3,...
             'densitySet',2,...
                      'T',600,...
                     'BW',10,...
             'enrichment',[0.032 0.032 0.032],...
                 'isBare',0,...
            'makeLibrary',1,...
         'printDiffParam',0,...
                     'gp',300,...
                    'BCs',[],...
              'makePlots',0,...
             'plotChoice',[1 1]);
%              'plotName',{});
%         p = struct('nDim',2,...
%                 'nGroups',3,...
%              'densitySet',2,...
%                       'T',600,...
%                      'BW',2.0246,...
%              'enrichment',[0.032 0.032 0.032],...
%                  'isBare',0,...
%             'makeLibrary',1,...
%          'printDiffParam',0,...
%                      'gp',300,...
%               'mapChoice',0,...
%                     'BCs',[],...
%               'makePlots',1,...
%              'plotChoice',[1 1]);
    end
    % DEFINE GEOMETRY
    g = makeGeometry(p);
end
% MAKE LIBRARY
if p.makeLibrary
    disp('Making library')
    L = makeLibrary(p);
    save XSLibrary.mat L
else
    if ~exist('XSLibrary.mat','var')
        disp('Loading library')
        load XSLibrary.mat
    end
end


% INTERPOLATE LIBRARY
disp('Interpolating library')
Li = interpolateLibrary(L,p);

% COLLAPSE LIBRARY
disp('Collapsing library')
Region = collapseLibrary(Li,g,p);

% ASSEMBLE DIFFUSION PARAMETERS
disp('Assembling diffusion parameters')
DiffParam = makeDiffParam(L,Region,g,p);

A = [ DiffParam{1}.xssink DiffParam{1}.vfission DiffParam{1}.diffco];

dlmwrite('ds2bw10.txt',A,'&');
    
end


%{
% mod for bare reactor
if p.isBare
    DiffParam{regrefl}.scatker = zeros(p.nGroups);
    DiffParam{regrefl}.xssink = zeros(1,p.nGroups);
    DiffParam{regrefl}.vfission = zeros(1,p.nGroups);
    DiffParam{regrefl}.diffco = 100000*ones(1,p.nGroups);
end

% mod for 1-group calculation
if p.nGroups == 1
    for regidx = 1:g.nRegions
        DiffParam{regidx}.xssink = Region(regidx).NBins(L.Bin(p.nGroups)).xss(L.MT(18)).xs...
                                 + Region(regidx).NBins(L.Bin(p.nGroups)).xss(L.MT(102)).xs;
        DiffParam{regidx}.vfission = 2.5*DiffParam{regidx}.vfission;
    end
end

% RUN DIFFUSION
if p.nDim == 1
    disp('Running 1D diffusion')
%     Results = struct('R',R,...
%                   'Flux',fluxOut,...
%                   'keff',keff,...
%              'powerFrac',powerFrac,...
%                'runtime',runtime,...
%            'nIterations',counter-1,...
%                  'Error',ErrrVector);
    Results = Diffusion1D(DiffParam,g,p);

elseif p.nDim == 2
    disp('Running 2D diffusion')
% Results = struct('Xmesh',Xmesh,...
%                  'Ymesh',Ymesh,...
%                   'Flux',fluxOut,...
%                   'keff',keff,...
%              'powerGrid',powerEff,...
%                'runtime',runtime,...
%            'nIterations',counter-1,...
%                  'Error',ErrrVector);
%              'powerFrac',powerFrac,...
%               'gpBuffer',gpBuffer

    if p.mapChoice == 0
        g2D = createMapSimple(g,p);
    elseif p.mapChoice == 4
        g2D = createMap4(g,p);
    elseif p.mapChoice == 5 % NOT READY YET
        g2D = createMap5(g,p);
    end
    Results = Diffusion2D(DiffParam,g2D,g,p);

end
disp('Done with diffusion!')

% PLOT RESULTS
if p.makePlots
    disp('Making plots')
    if p.nDim == 1
        
        if p.nGroups == 3
            Diffusion1D3gPlots(Results,g,p);
        elseif p.nGroups == 1
            Diffusion1D1gPlots(Results,g,p);
        end
        
    elseif p.nDim == 2
        
        if p.nGroups == 3
            
            Diffusion2D3gPlots(Results,g,p);
            
            if p.mapChoice == 4
                Diffusion2D3gA4Plots(Results,g,p);
            elseif p.mapChoice == 5 % NOT READY YET
%                 Diffusion2D3 
            end
            
        elseif p.nGroups == 1
            Diffusion2D1gPlots(Results,g,p);
        end
        
    end  
end
%}