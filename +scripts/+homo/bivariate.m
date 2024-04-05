function bivariate(savename, reactortype, homoset, mcnpxfluxM)
threshfacs = [.3:.2:5];
slopefacs = [.5:.2:4];
nT = length(threshfacs);
nS = length(slopefacs);
kcompare = zeros(nT, nS);
speccompare = cell(nT, nS);
specerror = zeros(nT, nS);
homocompare = cell(nT, nS);
for iT = 1:length(threshfacs)
    for iS = 1:length(slopefacs)
        [p, g] = scripts.homo.inputdef(reactortype, homoset, ...
            threshfacs(iT), slopefacs(iS));
        p.homogenizeThreshFac = threshfacs(iT); %threshfac;
        p.homogenizeSlopeFac = slopefacs(iS); %slopefac;

        [RM, pM, gM, LM] = vbudsii.Vbudsii(p, g);

        vbudsiifluxM = RM.Region(1).spectralFlux;
        vbudsiifluxM = vbudsiifluxM / sum(sum(vbudsiifluxM)) * 100;
        vbudsiikM = RM.Region(1).kInf;

        kcompare(iT,iS) = vbudsiikM;
        speccompare{iT,iS} = vbudsiifluxM;
        specerror(iT, iS) = sum(sum(abs(vbudsiifluxM/sum(vbudsiifluxM(:,1)) -...
            mcnpxfluxM/sum(mcnpxfluxM(:,1)))));
        homocompare{iT,iS} = RM.Region(1).Homo;
    end
end

save(savename, 'threshfacs', 'slopefacs', ...
    'kcompare', 'speccompare', 'specerror', 'homocompare')


