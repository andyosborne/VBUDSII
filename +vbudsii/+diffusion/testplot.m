% testplot

p1 = struct('nDim',1,...
    'nGroups',3,...
    'densitySet',2,...
    'T',600,...
    'BW',2.0246,...
    'enrichment',[0.032 0.032 0.032],...
    'isBare',0,...
    'makeLibrary',1,...
    'printDiffParam',0,...
    'gp',30,...
    'mapChoice',0,...
    'makePlots',1,...
    'plotChoice',[1 0 1],...
    'plotSave',[0 0 0],...
    'plotName',[figdir 'pg_gridspacing_0510']);

Resultz = main(p1);