
    import scripts.homo.*;
    import plotmachines.*;

    runs = plotmachines.FourCells('multicell4cellfastgeoff');

    nohomo = 'Method unmodified.';

    fastgeoffdescript = ['FAST 2. Pin 1: natural UO2, ' ...
        'Pin 2: PuO2. '];


    [pFG, gFG] = inputdef('fastgeoff', 'off');
    runs.AddRun('fastgeoff', pFG, gFG, [fastgeoffdescript nohomo]);
    runs.Run('fastgeoff', true);
    theLib = runs.LoadLib('fastgeoff');
    [mcnpxfluxFG, mcnpxkFG, mcnpxRFG] = ohAlsoMCNPX(...
        pFG, gFG, theLib, 'FG', true);
    runs.SetMcnpx('fastgeoff', mcnpxfluxFG, mcnpxkFG);
    runs.Load('fastgeoff');
    input2tablemulti(gFG, 'multifastgeoffinput');
    RFG = runs.Results('fastgeoff');
    vbudsiifluxFG = RFG.Region(1).spectralFlux * 100 / ...
        sum(RFG.Region(1).spectralFlux(:,1));
    mcnpxfluxFG = mcnpxfluxFG / sum(mcnpxfluxFG(:,1)) * 100;
    lowbound = 60;
    plotmachines.fluxandxsgeoff('flux', theLib.groupDef(lowbound:end), ...
        mcnpxfluxFG(lowbound:end,1), vbudsiifluxFG(lowbound:end,1), ...
        'printit', 'multifastgeoffuo2flux');
    plotmachines.fluxandxsgeoff('flux', theLib.groupDef(lowbound:end), ...
        mcnpxfluxFG(lowbound:end,3), vbudsiifluxFG(lowbound:end,3), ...
        'printit', 'multifastgeoffplutflux');

    theLib = runs.LoadLib('fastgeoff');

    disp(sprintf('our keff: %f, mcnpx keff: %f', RFG.Region(1).kInf,...
    mcnpxkFG));

    % Reaction rate plots.
    plotmachines.reactionrates.makeRRtables('multi', ...
        'MCNPX', RFG, mcnpxRFG, pFG, gFG, theLib);

    runs.printreport()

    VbudsiiVsMcnpx('multifastgeoffv2m', 'multifastgeoffv2m', ...
        gFG, theLib, RFG, mcnpxRFG, ...
        sum(sum(vbudsiifluxFG))/sum(sum(mcnpxfluxFG)));
