function g = initgeom()

g = struct('nRegions', [],...
          'regionDef', struct('name','',...
                            'nCells',    [],...
                        'relVolumes',[],...
                     'isFissionable',[],...
                           'cellDef',struct('name','',...
                                   'isFissionable',[],...
                                       'initZAIDs',[],...
                                'initNumDensities',[],...
                                     'initDensity',[],...
                                        'initTemp',[]...
    )));

end
