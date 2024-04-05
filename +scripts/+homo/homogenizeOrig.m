function homogenizeOrig
    set(0, 'DefaultAxesFontSize', 8);
    set(0, 'DefaultTextFontSize', 8);
    set(0, 'DefaultTextFontSize', 8);

    runvbudsii = false;

    % 1-make sure input is working
    % 1.5- make the new library
    % 2-tie to mcnpx runs
    % 3-create methods for making a report
    % 4-finally, tweak parameters, etc.
    % Create a nice report
    % Want a nicely-formatted four-pane figure, like what we had before, but
    % cleaner!
    % Show the same-fuel case first
    % Then show low-enriched uranium
    % In all cases make sure all four plots have the same axis scaling, etc.
    % - captions
    % - keff
    % - homogeneous thing
    %
    % plot error, redo runs for min error, send to mark.
    % Possible bug: why thermal hump gets better when moving to individual
    % homogeneity.
    import scripts.homo.*;
    import plotmachines.*;

    runs = plotmachines.FourCells('homogenizestudy');

    %
    [pE, gE] = inputdef('equi', 'off');
    runs.AddRun('equi', pE, gE, ...
        'Both fuel pins are the same enriched UOX. Method unmodified.');
    if runvbudsii
        runs.Run('equi', true);
    end
    theLib = runs.LoadLib('equi');
    [mcnpxfluxE, mcnpxkE, mcnpxRE] = ohAlsoMCNPX(...
        pE, gE, theLib, 'E', false);
    runs.SetMcnpx('equi', mcnpxfluxE, mcnpxkE);
    runs.Load('equi');

    nohomo = 'Method unmodified. ';
    orighomo = 'Homogeneity: original scheme. ';
    indivhomo = 'Homogeneity: 2nd ''indiv'' scheme. ';

    puredescript = ['THERMAL 1. '...
        'Pin 1: 3\\%% enriched UOX, Pin 2: 0\\%% enriched UOX. '];
    imfdescript = ['THERMAL 2. ' ...
        'Pin 1: 3\\%% enriched UOX, Pin 2: IMF. '];
    fastdescript = ['FAST 1. ' ...
        'Pin 1: 24\\%% enriched UOX, Pin 2: 12\\%% enriched UOX'];
    pufastdescript = ['FAST 2. Pin 1: 24\\%% enriched UOX, ' ...
        'Pin 2: Pu(239)O2. '];


    % Pure
    % ----
    [pP, gP] = inputdef('pure', 'off');
    runs.AddRun('pure', pP, gP, [puredescript nohomo]);
    if runvbudsii
        runs.Run('pure', false);
    end
    [mcnpxfluxP, mcnpxkP, mcnpxRP] = ohAlsoMCNPX(...
        pP, gP, theLib, 'P', false);
    runs.SetMcnpx('pure', mcnpxfluxP, mcnpxkP);
    runs.Load('pure');


    [pPHO, gPHO] = inputdef('pure', 'orig');
    runs.AddRun('pureorig', pPHO, gPHO, [puredescript orighomo]);
    if runvbudsii
        runs.Run('pureorig', false);
    end
    [mcnpxfluxPHO, mcnpxkPHO, mcnpxRPHO] = ohAlsoMCNPX(...
        pPHO, gPHO, theLib, 'PHO', false);
    runs.SetMcnpx('pureorig', mcnpxfluxPHO, mcnpxkPHO);
    commentPHO = ['In this method, a single homogenization factor ' ...
        'is used for all cells. ' ...
        'The multiplication factor is slightly better ' ...
        'using this simple homogeneity method, but the thermal ' ...
        'peak in the pure UO2 moderator increases, whereas it had been ' ...
        'correct before using homogenization.'];
    runs.Load('pureorig', commentPHO);

    [pPHI, gPHI] = inputdef('pure', 'indiv');
    pPFHI.homogenizeThreshFac = .5; %threshfac;
    pPFHI.homogenizeSlopeFac = 1.1; %slopefac;
    %The key improvement we see here
    runs.AddRun('pureindiv', pPHI, gPHI, [puredescript indivhomo]);
    if runvbudsii
        runs.Run('pureindiv', false);
    end
    [mcnpxfluxPHI, mcnpxkPHI, mcnpxRPHI] = ohAlsoMCNPX(...
        pPHI, gPHI, theLib, 'PHO', false);
    runs.SetMcnpx('pureindiv', mcnpxfluxPHI, mcnpxkPHI);
    commentPHI = ['In this case, each cell has its own homogenization ' ...
        'factor, and the resulting PI is obtained in a mannner ' ...
        'different than that used for Figure 3. The multiplication ' ...
        'factor is worse than for the simple homogeneity method, ' ...
        'but the flux in the pure UOX thermal now again matches the ' ...
        'MCNPX result (as it did when not doing homogenization).' ...
        ' Note that this output was generated using the two parameters ' ...
        'for homogeneity that yielded the least error.'];
    runs.Load('pureindiv', commentPHI);

    % IMF
    % ---
    [pI, gI] = inputdef('imf', 'off');
    runs.AddRun('imf', pI, gI, [imfdescript nohomo]);
    if runvbudsii
        runs.Run('imf', false);
    end
    [mcnpxfluxI, mcnpxkI, mcnpxRI] = ohAlsoMCNPX(...
        pI, gI, theLib, 'I', false);
    runs.SetMcnpx('imf', mcnpxfluxI, mcnpxkI);
    runs.Load('imf');
    input2tablemulti(gI, 'multiimfinput');
    RI = runs.Results('imf');
    vbudsiifluxI = RI.Region(1).spectralFlux * 100 / ...
        sum(RI.Region(1).spectralFlux(:,1));
    mcnpxfluxI = mcnpxfluxI / sum(mcnpxfluxI(:,1)) * 100;
    plotmachines.fluxandxsgeoff('flux', theLib.groupDef, ...
        mcnpxfluxI(:,1), vbudsiifluxI(:,1), ...
        'printit', 'multiimfuo2flux');
    plotmachines.fluxandxsgeoff('flux', theLib.groupDef, ...
        mcnpxfluxI(:,3), vbudsiifluxI(:,3), ...
        'printit', 'multiimfimfflux');


    [pIHO, gIHO] = inputdef('imf', 'orig');
    runs.AddRun('imforig', pIHO, gIHO, [imfdescript orighomo]);
    if runvbudsii
        runs.Run('imforig', false);
    end
    [mcnpxfluxIHO, mcnpxkIHO, mcnpxRIHO] = ohAlsoMCNPX(...
        pIHO, gIHO, theLib, 'IHO', false);
    runs.SetMcnpx('imforig', mcnpxfluxIHO, mcnpxkIHO);
    commentIHO = ['Again, one homogenization factor is used for all cells.' ...
        ' The multiplication factor improves, but the flux does not change.'];
    runs.Load('imforig', commentIHO);
    RIHO = runs.Results('imforig');
    vbudsiifluxIHO = RIHO.Region(1).spectralFlux * 100 / ...
        sum(RIHO.Region(1).spectralFlux(:,1));
    mcnpxfluxIHO = mcnpxfluxIHO / sum(mcnpxfluxIHO(:,1)) * 100;
    plotmachines.fluxandxsgeoff('flux', theLib.groupDef, ...
        mcnpxfluxIHO(:,1), vbudsiifluxIHO(:,1), ...
        'printit', 'multiimfhomouo2flux');
    plotmachines.fluxandxsgeoff('flux', theLib.groupDef, ...
        mcnpxfluxIHO(:,3), vbudsiifluxIHO(:,3), ...
        'printit', 'multiimfhomoimfflux');
    
    [pIHI, gIHI] = inputdef('imf', 'indiv');
    pPFHI.homogenizeThreshFac = .3; %threshfac;
    pPFHI.homogenizeSlopeFac = .5; %slopefac;
    runs.AddRun('imfindiv', pIHI, gIHI, [imfdescript indivhomo]);
    if runvbudsii
        runs.Run('imfindiv', false);
    end
    [mcnpxfluxIHI, mcnpxkIHI, mcnpxRIHI] = ohAlsoMCNPX(...
        pIHI, gIHI, theLib, 'IHO', false);
    runs.SetMcnpx('imfindiv', mcnpxfluxIHI, mcnpxkIHI);
    commentIHI = ['Multiplication factor gets worse and the flux does ' ...
        'not improve.' ...
        ' Note that this output was generated using the two parameters ' ...
        'for homogeneity that yielded the least error.'];
    runs.Load('imfindiv', commentIHI);

    % Fast
    % ----
    [pF, gF] = inputdef('fast', 'off');
    runs.AddRun('fast', pF, gF, [fastdescript nohomo]);
    if runvbudsii
        runs.Run('fast', false);
    end
    [mcnpxfluxF, mcnpxkF, mcnpxRF] = ohAlsoMCNPX(...
        pF, gF, theLib, 'F', false);
    runs.SetMcnpx('fast', mcnpxfluxF, mcnpxkF);
    runs.Load('fast');

    [pFHO, gFHO] = inputdef('fast', 'orig');
    runs.AddRun('fastorig', pFHO, gFHO, [fastdescript orighomo]);
    if runvbudsii
        runs.Run('fastorig', false);
    end
    [mcnpxfluxFHO, mcnpxkFHO, mcnpxRFHO] = ohAlsoMCNPX(...
        pFHO, gFHO, theLib, 'FHO', false);
    runs.SetMcnpx('fastorig', mcnpxfluxFHO, mcnpxkFHO);
    commentFHO = ['As expected, homogeneity has no appreciable effect on ' ...
        'fast reactor results.'];
    runs.Load('fastorig', commentFHO);

    [pFHI, gFHI] = inputdef('fast', 'indiv');
    runs.AddRun('fastindiv', pFHI, gFHI, [fastdescript indivhomo]);
    if runvbudsii
        runs.Run('fastindiv', false);
    end
    [mcnpxfluxFHI, mcnpxkFHI, mcnpxRFHI] = ohAlsoMCNPX(...
        pFHI, gFHI, theLib, 'FHO', false);
    runs.SetMcnpx('fastindiv', mcnpxfluxFHI, mcnpxkFHI);
    runs.Load('fastindiv');

    % Fast plut.
    % ----------
    [pPF, gPF] = inputdef('fastplut', 'off');
    runs.AddRun('pufast', pPF, gPF, [pufastdescript nohomo]);
    if runvbudsii
        runs.Run('pufast', false);
    end
    [mcnpxfluxPF, mcnpxkPF, mcnpxRPF] = ohAlsoMCNPX(...
        pPF, gPF, theLib, 'PF', false);
    runs.SetMcnpx('pufast', mcnpxfluxPF, mcnpxkPF);
    runs.Load('pufast');
    input2tablemulti(gPF, 'multifastplutinput');
    RPF = runs.Results('pufast');
    vbudsiifluxPF = RPF.Region(1).spectralFlux * 100 / ...
        sum(RPF.Region(1).spectralFlux(:,1));
    mcnpxfluxPF = mcnpxfluxPF / sum(mcnpxfluxPF(:,1)) * 100;
    lowbound = 60;
    plotmachines.fluxandxsgeoff('flux', theLib.groupDef(lowbound:end), ...
        mcnpxfluxPF(lowbound:end,1), vbudsiifluxPF(lowbound:end,1), ...
        'printit', 'multifastplutuo2flux');
    plotmachines.fluxandxsgeoff('flux', theLib.groupDef(lowbound:end), ...
        mcnpxfluxPF(lowbound:end,3), vbudsiifluxPF(lowbound:end,3), ...
        'printit', 'multifastplutplutflux');

    [pPFHO, gPFHO] = inputdef('fastplut', 'orig');
    runs.AddRun('pufastorig', pPFHO, gPFHO, [pufastdescript orighomo]);
    if runvbudsii
        runs.Run('pufastorig', false);
    end
    [mcnpxfluxPFHO, mcnpxkPFHO, mcnpxRPFHO] = ohAlsoMCNPX(...
        pPFHO, gPFHO, theLib, 'PFHO', false);
    runs.SetMcnpx('pufastorig', mcnpxfluxPFHO, mcnpxkPFHO);
    runs.Load('pufastorig');

    [pPFHI, gPFHI] = inputdef('fastplut', 'indiv');
    runs.AddRun('pufastindiv', pPFHI, gPFHI, [pufastdescript indivhomo]);
    if runvbudsii
        runs.Run('pufastindiv', false);
    end
    [mcnpxfluxPFHI, mcnpxkPFHI, mcnpxRPFHI] = ohAlsoMCNPX(...
        pPFHI, gPFHI, theLib, 'PFHO', false);
    runs.SetMcnpx('pufastindiv', mcnpxfluxPFHI, mcnpxkPFHI);
    runs.Load('pufastindiv');

    %bivariate('bivariatePHI', 'pure', 'indiv', mcnpxfluxPHI);
    %bivariate('bivariateFHI', 'fast', 'indiv', mcnpxfluxFHI);
    %bivariate('bivariatePFHI', 'fastplut', 'indiv', mcnpxfluxPFHI);
    %bivariate('bivariateIHI', 'imf', 'indiv', mcnpxfluxIHI);
        
    save('mcnpxkeffs', 'mcnpxkP', 'mcnpxkI', 'mcnpxkF', 'mcnpxkPF');
    runs.printreport()

end



    %[RP, pP, gP, LP] = vbudsii.Vbudsii(pP, gP);
    %[RE, pE, gE, LE] = vbudsii.Vbudsii(pE, gE);
    %[RPHO, pPHO, gPHO, LPHO] = vbudsii.Vbudsii(pPHO, gPHO);
    %[RF, pF, gF, LF] = vbudsii.Vbudsii(pF, gF);
    %[RFHO, pFHO, gFHO, LFHO] = vbudsii.Vbudsii(pFHO, gFHO);
    %[RI, pI, gI, LI] = vbudsii.Vbudsii(pI, gI);
    %[RIHO, pIHO, gIHO, LIHO] = vbudsii.Vbudsii(pIHO, gIHO);




