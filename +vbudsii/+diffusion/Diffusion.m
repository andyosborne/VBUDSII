function Region = Diffusion(L, p, g, Region)
%MAIN Run difussion with a parameter structure p and a geometry structure g.
% Description: Makes a Library structure, interpolates it for a temperature
% T and Bonderenko self-shielding factor BW, collapses library according to
% geometry g, creates diffusion parameters, and outputs Results structure
% from diffusion.
%
% USE:  Results = MAIN()
%       Results = MAIN(p)
%       Results = MAIN(p,g)
%
% NOTES: The typical use of this function is to run it in a script that
% defines a parameter structure p, calls MAIN(p), and operates on the
% Results in the way that pleases it.
%
% EXAMPLES:
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110517 cld72   cleaned up and formatted
%
% FUTURE UPDATES:
%   1- work without defining its own library; integrate into remainder of
%   code
%
% DEPENDENCIES:
%   makeLibrary
%   makeGeometry
%   interpolateLibrary
%   collapseLibrary
%   makeDiffParam
%   Diffusion1D
%   Diffusion2D
%   Diffusion1D1gPlots
%   Diffusion1D3gPlots
%   Diffusion2D1gPlots
%   Diffusion2D3gPlots
%   Diffusion2D3gA4Plots
%

if p.verbose
    disp('Entering Diffusion');
end

%{
fprintf('\n!!!New Run!!! at %s\n',datestr(now,31))
close all;

% manage variable input
if length(varargin) == 2
    % both g and p are inputs
    g = varargin{1};
    p = varargin{2};
elseif length(varargin) == 1 || isempty(varargin)
    % only p is input, or neither p or g is input.
    if length(varargin) == 1
        % only p is input
        p = varargin{1};
    elseif isempty(varargin)
        % DEFINE PARAMETERS
        p = struct('nDim',1,...
                'nGroups',1,...
             'densitySet',2,...
                      'T',600,...
                     'BW',2.0246,...
             'enrichment',[0.032 0.032 0.032],...
                 'isBare',1,...
            'makeLibrary',0,...
         'printDiffParam',0,...
                     'gp',300,...
                    'BCs',[],...
              'makePlots',1,...
             'plotChoice',[1 1 1],...
               'plotSave',[0 0 0],...
               'plotName','test');
    end
    % DEFINE GEOMETRY
    g = makeGeometry(p);
end
% MAKE LIBRARY
if p.makeLibrary
    % parameter
    disp('Making library')
    L = makeLibrary(p);
    save XSLibrary.mat L
else
    if ~exist('XSLibrary.mat','var')
        disp('Loading library')
        load XSLibrary.mat
    end
end

% ONLY 3-GROUP DIFFUSION IS SUPPORTED RIGHT NOW.
% INTERPOLATE LIBRARY
disp('Interpolating library')
Li = interpolateLibrary(L,p);
%}

%% ERROR CHECK
if ~(p.nFewGroups == 1 || p.nFewGroups == 3)
    % The few number of groups is not 1 or 3. Diffusion can only work with 1 or
    % 3 few groups.
    error(['ERROR: The value of p.nFewGroups is ' num2str(p.nFewGroups) ...
           ', but Diffusion currently only supports either 1 or 3' ...
           ' few groups']);
end

% Collapse fine group cell cross sections to few group region cross sections.
Region = ResolveRegionXS(L, p, g, Region);

% Assemble diffusion parameters.
DiffParam = MakeDiffParam(L, p, g, Region);

% modify diffparams for bare reactor
if p.isBare
    regrefl = 4; % region index for reflector region
    DiffParam{regrefl}.scatker = zeros(p.nGroups);
    DiffParam{regrefl}.xssink = zeros(1,p.nGroups);
    DiffParam{regrefl}.vfission = zeros(1,p.nGroups);
    DiffParam{regrefl}.diffco = 100000*ones(1,p.nGroups);
end

% modify diffparams for 1-group calculation
if p.nGroups == 1
    for regidx = 1:g.nRegions
        DiffParam{regidx}.xssink = ...
            Region(regidx).few(L.MT(18)).value...
               + Region(regidx).few(L.MT(102)).value;
        DiffParam{regidx}.vfission = 2.5*DiffParam{regidx}.vfission;
    end
end

% RUN DIFFUSION
if p.nDim == 1
%     Results = struct('R',R,...
%                   'Flux',fluxOut,...
%                   'keff',keff,...
%              'powerFrac',powerFrac,...
%                'runtime',runtime,...
%            'nIterations',counter-1,...
%                  'Error',ErrrVector);
    Results = Diffusion1D(DiffParam,g,p);

elseif p.nDim == 2
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
        % g2D = createMap5(g,p);
    end
    Results = Diffusion2D(DiffParam,g2D,g,p);

end

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
                % Diffusion2D3gA5Plots(Results,g,p);
            end
        elseif p.nGroups == 1
            Diffusion2D1gPlots(Results,g,p);
        end
    end
end

if p.verbose
    disp('Leaving Diffusion');
end
