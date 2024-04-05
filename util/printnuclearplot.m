function g = printnuclearplot(fname, datatype, legends, groupdefs, xss, ...
    varargin)

    ip = inputParser;

    validDatatype = {'micro', 'macro', 'flux', 'errormicro', 'errormacro', ...
        'errorflux', 'pi'};
    checkDatatype = @(x) any(validatestring(x, validDatatype));

    defaultFdir = '';

    defaultEvince = false;

    defaultBinit = false;

    defaultMytitle = '';

    defaultXbelow = '';

    defaultYlims = [];

    defaultGroupdefs2 = {};

    defaultXss2 = {};

    defaultYlog = false;

    defaultS0bounds = false;

    defaultSmall = false;

    defaultLegendOn = true;

    addRequired(ip, 'fname', @ischar);
    addRequired(ip, 'datatype', checkDatatype);
    addRequired(ip, 'legends', @iscell);
    addRequired(ip, 'groupdefs');
    addRequired(ip, 'xss', @iscell);

    addOptional(ip, 'fdir', defaultFdir, @isdir);
    addOptional(ip, 'evince', defaultEvince, @islogical);
    addOptional(ip, 'binit', defaultBinit, @islogical);
    addOptional(ip, 'mytitle', defaultMytitle, @ischar);
    addOptional(ip, 'xbelow', defaultXbelow, @ischar);
    addOptional(ip, 'ylims', defaultYlims, @isnumeric);
    addOptional(ip, 'ylog', defaultYlog, @islogical);
    addOptional(ip, 'groupdefs2', defaultGroupdefs2);
    addOptional(ip, 'xss2', defaultXss2, @iscell);
    addOptional(ip, 's0bounds', defaultS0bounds);
    addOptional(ip, 'small', defaultSmall, @islogical);
    addOptional(ip, 'legendOn', defaultLegendOn, @islogical);

    parse(ip, fname, datatype, legends, groupdefs, xss, varargin{:});

    fdir = ip.Results.fdir;
    evince = ip.Results.evince;
    binit = ip.Results.binit;
    mytitle = ip.Results.mytitle;
    xbelow = ip.Results.xbelow;
    ylims = ip.Results.ylims;
    ylog = ip.Results.ylog;
    groupdefs2 = ip.Results.groupdefs2;
    xss2 = ip.Results.xss2;
    s0bounds = ip.Results.s0bounds;
    small = ip.Results.small;
    legendOn = ip.Results.legendOn;

if xor(isempty(groupdefs2), isempty(xss2))
    error(['Either groupdefs2 and xss2 must both be defined, or neither ' ...
        'are defined.'])
end

if ~isempty(fdir) && ~strcmp(fdir(end), filesep)
    fdir = [fdir filesep];
end
fdirname = [fdir fname];

nplots = length(xss);
lengthX = zeros(nplots,1);
if iscell(groupdefs)
    for i = 1:nplots
        lengthX(i) = length(groupdefs{i}) - 1;
        if binit == 1
            XaxisXS{i} = zeros(2*lengthX(i),1);
            XaxisXS{i} = groupdefs{i}(1);
            for j = 2:lengthX(i)

                XaxisXS{i}(2*j-2:2*j-1) = groupdefs{i}(j)*[1; 1];
            end
            XaxisXS{i}(2*lengthX(i)) = groupdefs{i}(end);
        else
            XaxisXS{i} = groupdefs{i}(1:end-1)';
        end
    end
else
    for i = 1:nplots
        lengthX(i) = length(groupdefs) - 1;
        if binit == 1
            XaxisXS{i} = zeros(2*lengthX(i),1);
            XaxisXS{i} = groupdefs(1);
            for j = 2:lengthX(i)

                XaxisXS{i}(2*j-2:2*j-1) = groupdefs(j)*[1; 1];
            end
            XaxisXS{i}(2*lengthX(i)) = groupdefs(end);
        else
            XaxisXS{i} = groupdefs(1:end-1)';
        end
    end
end
% prepare for bin-wise plotting.
plotn = cell(nplots,1);
minxval = inf;
maxxval = -inf;
minyval = inf;
maxyval = -inf;
for i = 1:nplots
    if binit == 1
        plotn{i} = zeros(lengthX(i)*2,1);
        plotn{i}(1:2:2*lengthX(i)-1) = xss{i};
        plotn{i}(2:2:2*lengthX(i)) = xss{i};
    else
        plotn{i} = xss{i};
    end

    maxx = max(XaxisXS{i});
    if maxx(1) > maxxval
        maxxval = maxx(1);
    end
    minx = min(XaxisXS{i});
    if minx(1) < minxval
        minxval = minx(1);
    end
    maxy = max(plotn{i});
    if maxy(1) > maxyval
        maxyval = maxy(1);
    end
    miny = min(plotn{i});
    if miny(1) < minyval
        minyval = miny(1);
    end
end

if ~isempty(groupdefs2)
    nplots2 = length(xss2);
    lengthX2 = zeros(nplots2,1);
    if iscell(groupdefs2)
        for i = 1:nplots2
            lengthX2(i) = length(groupdefs2{i}) - 1;
            if binit == 1
                XaxisXS2{i} = zeros(2*lengthX2(i),1);
                XaxisXS2{i} = groupdefs2{i}(1);
                for j = 2:lengthX2(i)

                    XaxisXS2{i}(2*j-2:2*j-1) = groupdefs2{i}(j)*[1; 1];
                end
                XaxisXS2{i}(2*lengthX2(i)) = groupdefs2{i}(end);
            else
                XaxisXS2{i} = groupdefs2{i}(1:end-1)';
            end
        end
    else
        for i = 1:nplots2
            lengthX2(i) = length(groupdefs2) - 1;
            if binit == 1
                XaxisXS2{i} = zeros(2*lengthX2(i),1);
                XaxisXS2{i} = groupdefs2(1);
                for j = 2:lengthX2(i)

                    XaxisXS2{i}(2*j-2:2*j-1) = groupdefs2(j)*[1; 1];
                end
                XaxisXS2{i}(2*lengthX2(i)) = groupdefs2(end);
            else
                XaxisXS2{i} = groupdefs2(1:end-1)';
            end
        end
    end
    % prepare for bin-wise plotting.
    plotn2 = cell(nplots2,1);
    for i = 1:nplots2
        if binit == 1
            plotn2{i} = zeros(lengthX2(i)*2,1);
            plotn2{i}(1:2:2*lengthX2(i)-1) = xss2{i};
            plotn2{i}(2:2:2*lengthX2(i)) = xss2{i};
        else
            plotn2{i} = xss2{i};
        end

        if max(XaxisXS2{i}) > maxxval
            maxxval = max(XaxisXS2{i});
        end
        if min(XaxisXS2{i}) < minxval
            minxval = min(XaxisXS2{i});
        end
        if max(plotn2{i}) > maxyval
            maxyval = max(plotn2{i});
        end
        if min(plotn2{i}) < minyval
            minyval = min(plotn2{i});
        end
    end
end

fid = fopen([fdirname '.tex'], 'w');

fprintf(fid,'\\documentclass{minimal}\n');
fprintf(fid,'\\usepackage{silence}\n');
fprintf(fid,'\\usepackage{times}\n');
fprintf(fid,'\\usepackage{helvet}\n');
fprintf(fid,'\\renewcommand{\\familydefault}{\\sfdefault}\n');
fprintf(fid,'\\usepackage{sfmath}\n');
fprintf(fid,'\\usepackage{pgfplots}\n');
fprintf(fid,'\\usepackage[active,tightpage]{preview}\n');
fprintf(fid,'\\PreviewEnvironment{tikzpicture}\n');
fprintf(fid,'\\setlength\\PreviewBorder{5pt}\n');
fprintf(fid,'\\begin{document}\n');
fprintf(fid,'\\WarningsOff\n');
fprintf(fid,'\\begin{tikzpicture}\n');
if ylog
    fprintf(fid,'\\begin{loglogaxis}[\n');
else
    fprintf(fid,'\\begin{semilogxaxis}[\n');
end
%fprintf(fid,'   y=3,\n');
%fprintf(fid,'   x=5,\n');
if legendOn
%fprintf(fid,'   legend style={at={(.05,.95)},anchor=north west},\n');
%fprintf(fid,'   legend style={at={(.05,.05)},anchor=south west},\n');
%fprintf(fid,'   legend style={at={(.05,.05)},anchor=south west},\n');
fprintf(fid,'   legend cell align=left,\n');
if strcmp(datatype, 'flux')
    fprintf(fid,'   legend pos={north west},\n');
else
    fprintf(fid,'   legend pos={south west},\n');
end
end
fprintf(fid,'   cycle list name=exotic,\n');
if small
    fprintf(fid,'   width=2.8in,\n');
    fprintf(fid,'   height=2.5in,\n');
else
    fprintf(fid,'   width=6in,\n');
    fprintf(fid,'   height=2.5in,\n');
end
fprintf(fid,'   xmin=%.4e,\n',minxval);
fprintf(fid,'   xmax=%.4e,\n',maxxval);
if ~isempty(ylims)
    fprintf(fid,'   ymin=%.4e,\n',ylims(1));
    fprintf(fid,'   ymax=%.4e,\n',ylims(2));
else
    ylims(1) = minyval;
    ylims(2) = maxyval;
end
fprintf(fid,'   xlabel={Energy (eV) \\\\ %s},\n',xbelow);
fprintf(fid,'   xlabel style={text width=6cm, text centered},\n');
if strcmp(datatype, 'macro')
    myylabel = '$\Sigma$ (cm$^{-1}$)';
elseif strcmp(datatype, 'micro')
    myylabel = '$\sigma$ (b)';
elseif strcmp(datatype, 'flux')
    myylabel = '$\phi$ (n cm$^{-2}$ s$^{-1}$)';
elseif strcmp(datatype, 'errormacro')
    myylabel = 'error in $\Sigma$ (cm$^{-1}$)';
elseif strcmp(datatype, 'errormicro')
    myylabel = 'error in $\sigma$ (b)';
elseif strcmp(datatype, 'errorflux')
    myylabel = 'error in $\phi$ (n cm$^{-2}$ s$^{-1}$)';
elseif strcmp(datatype, 'pi')
    myylabel = '$\Pi(j<-i)$';
end
fprintf(fid,'   ylabel=%s,\n',myylabel);
fprintf(fid,'   title=%s,\n',mytitle);
fprintf(fid,'   no markers,\n');
fprintf(fid,']\n');
%col0rs = [0 0 255;
%          0 128 0;
%          255 0 0;
%          0 191 191;
%          191 0 191;
%          191 191 0;
%          64 64 64];
col0rs = [79 129 189;
          192 80 77;
          155 187 89;
          128 100 162;
          247 150 70;
          0 176 240;
          0 0 255;
         0 128 0;
         255 0 0;
         0 191 191;
         191 0 191;
         191 191 0;
         64 64 64];
if ~isempty(groupdefs2)
    for i = 1:nplots2
        fprintf(fid,'\\addplot[line width=1pt, color={rgb,255:red,%i;green,%i;blue,%i}] coordinates{\n',min(255, round(1.3*col0rs(i,:))));
        for j = 1:length(XaxisXS2{i})
            fprintf(fid,'   (%.4e,%.4e)\n',XaxisXS2{i}(j),plotn2{i}(j));
        end
        fprintf(fid,'};\n');
    end
end
for i = 1:nplots
    fprintf(fid,['\\addplot[line width=1pt, densely dotted,color={rgb,255:red,%i;green,%i;blue,%i}]'...
        'coordinates{\n'],col0rs(i,:));
    for j = 1:length(XaxisXS{i})
        fprintf(fid,'   (%.4e,%.4e)\n',XaxisXS{i}(j),plotn{i}(j));
    end
    fprintf(fid,'};\n');
end
fprintf(fid,'\\addplot[black!50,dashed] coordinates{\n');
fprintf(fid,'(1,%.4e)\n',min(minyval,ylims(1)));
fprintf(fid,'(1,%.4e)\n',max(maxyval,ylims(2)));
fprintf(fid,'};\n');
fprintf(fid,'\\addplot[black!50,dashed] coordinates{\n');
fprintf(fid,'(100e3,%.4e)\n',min(minyval,ylims(1)));
fprintf(fid,'(100e3,%.4e)\n',max(maxyval,ylims(2)));
fprintf(fid,'};\n');
if legendOn
fprintf(fid,'\\legend{');
if ~isempty(groupdefs2)
    pl0ts = [plotn2;plotn];
    legends = [legends(nplots+1:end) legends(1:nplots)];
else
    pl0ts = plotn;
end
for i = 1:length(legends)-1
    if any(pl0ts{i})
        fprintf(fid,'{%s},\n',legends{i});
    else
        fprintf(fid,'{%s},\n',[legends{i} '-0']);
        %warning('A data array is all zeros and is not being plotted.')
    end
end
if any(pl0ts{end})
    fprintf(fid,'{%s}',legends{end});
else
    fprintf(fid,'{%s}',[legends{end} '-0']);
    %warning('A data array is all zeros and is not being plotted.')
end
fprintf(fid,'}\n',legends{end});
end
if ylog
    fprintf(fid,'\\end{loglogaxis}\n');
else
    fprintf(fid,'\\end{semilogxaxis}\n');
end
fprintf(fid,'\\end{tikzpicture}\n');
fprintf(fid,'\\end{document}\n');
fclose(fid);

%if strcmp(fdir, '')
%    system(sprintf('pdflatex %s.tex', fdirname));
%else
%    system(sprintf('pdflatex -output-directory=%s %s.tex',fdir, fdirname));
%end
if strcmp(fdir, '')
    system(sprintf('rubber --pdf %s.tex', fdirname));
else
    system(sprintf('rubber --into=%s --pdf %s.tex',fdir, fdirname));
end
if evince
    system(sprintf('evince %s.pdf',fdirname));
end

% cleanup
delete([fdirname '.tex']);
delete([fdirname '.log']);
delete([fdirname '.aux']);



