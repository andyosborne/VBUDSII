function [mcnpxfluxE mcnpxkE mcnpxRE] = ohAlsoMCNPX(pE, gE, LE, name, dorun)
    fprintf('ohAlsoMCNPX the run %s.\n', name);
    mcnpxstr = ['mcnpx' name];
    mcnpxE = mcnpxwrap4.McnpxWrap4(mcnpxstr, ...
        fullfile('+homo4reactor', mcnpxstr), mcnpxstr);
    mcnpxfluxE = mcnpxE.RunAndGetFlux(pE, gE, dorun);
    mcnpxfluxE = mcnpxfluxE / sum(sum(mcnpxfluxE)) * 100;
    mcnpxkE = mcnpxE.GetMultiplicationFactor();
    mcnpxRE = mcnpxE.GetResults(LE);
end

