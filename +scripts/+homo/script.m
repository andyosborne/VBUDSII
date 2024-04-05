
matnames = {'PHI', 'IHI', 'FHI', 'PFHI'};

load('mcnpxkeffs');


mcnpxkeffs = [mcnpxkP mcnpxkI mcnpxkF mcnpxkPF];
if 0
    for idx = 1:2
        load(sprintf('bivariate%s', matnames{idx}));
        figure;
        theError = abs(kcompare - mcnpxkeffs(idx));
        surf(slopefacs, threshfacs, theError);

        [junk minidx] = min(theError(:));
        [nR nC] = size(theError);
        [iminx iminy] = ind2sub(size(theError), minidx);
        theOne = junk * ones(nR, nC);
        theOne(iminx, iminy) = theError(iminx, iminy) * 1.2;

        disp(sprintf('slope factor. min: %f', slopefacs(iminx)));
        disp(sprintf('thresh factor. min: %f', threshfacs(iminy)));
        xlabel(sprintf('slope factor. min: %f', slopefacs(iminx)));
        ylabel(sprintf('thresh factor. min: %f', threshfacs(iminy)));
        title(matnames{idx});
        zlabel('keff error abs(V - M)');
        disp(matnames{idx});

        hold on;
        surf(slopefacs, threshfacs, theOne);


        % 1 - find minimum error for pure and imf
        % 2 - run simulations at those values
        % 3 - analyze those simulations (in the morning)
        % 4 - get working fastgeoff results
        % 5 - reaction rates.
        % 6 - enroll in courses.

    end
end


for idx = 1:2
    load(sprintf('bivariate%s', matnames{idx}));
    figure;
    theError = specerror; %abs(kcompare - mcnpxkeffs(idx));
    surf(slopefacs, threshfacs, theError);

    [junk minidx] = min(theError(:));
    [nR nC] = size(theError);
    [iminx iminy] = ind2sub(size(theError), minidx);
    theOne = junk * ones(nR, nC);
    theOne(iminx, iminy) = theError(iminx, iminy) * 1.2;

    disp(sprintf('slope factor. min: %f', slopefacs(iminx)));
    disp(sprintf('thresh factor. min: %f', threshfacs(iminy)));
    xlabel(sprintf('slope factor. min: %f', slopefacs(iminx)));
    ylabel(sprintf('thresh factor. min: %f', threshfacs(iminy)));
    title(matnames{idx});
    zlabel('spectrum error abs(V - M)');
    disp(matnames{idx});

    %hold on;
    %surf(slopefacs, threshfacs, theOne);
end
