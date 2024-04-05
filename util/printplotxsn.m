function g = printplotxsn(fdir,fname,evince,binit,macro,ylims,mytitle,legends,xbelow,groupdefs,xss,varargin)
if length(varargin) == 2
    groupdefs2 = varargin{1};
    xss2 = varargin{2};
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

    if max(XaxisXS{i}) > maxxval
        maxxval = max(XaxisXS{i});
    end
    if min(XaxisXS{i}) < minxval
        minxval = min(XaxisXS{i});
    end
    if max(plotn{i}) > maxyval
        maxyval = max(plotn{i});
    end
    if min(plotn{i}) < minyval
        minyval = min(plotn{i});
    end
end

if length(varargin) == 2
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
fprintf(fid,'\\begin{semilogxaxis}[\n');
%fprintf(fid,'   y=3,\n');
%fprintf(fid,'   x=5,\n');
fprintf(fid,'   legend style={at={(.05,.95)},anchor=north west},\n');
%fprintf(fid,'   legend style={at={(.05,.05)},anchor=south west},\n');
fprintf(fid,'   cycle list name=exotic,\n');
fprintf(fid,'   width=6in,\n');
fprintf(fid,'   height=2.5in,\n');
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
if macro == 1
    myylabel = '$\Sigma$ (cm$^{-1}$)';
elseif macro == 0
    myylabel = '$\sigma$ (b)';
elseif macro == 2
    myylabel = '$\phi$ (n cm$^{-2}$ s$^{-1}$)';
end
fprintf(fid,'   ylabel=%s,\n',myylabel);
fprintf(fid,'   title=%s,\n',mytitle);
fprintf(fid,'   no markers,\n');
fprintf(fid,']\n');
col0rs = [0 0 255;
          0 128 0;
          255 0 0;
          0 191 191;
          191 0 191;
          191 191 0;
          64 64 64];
col0rs = [79 129 189;
          192 80 77;
          155 187 89;
          128 100 162;
          247 150 70;
          0 176 240];
if length(varargin) == 2
    for i = 1:nplots2
        fprintf(fid,'\\addplot[line width=1.5pt,color={rgb,255:red,%i;green,%i;blue,%i}] coordinates{\n',col0rs(i,:));
        for j = 1:length(XaxisXS2{i})
            fprintf(fid,'   (%.4e,%.4e)\n',XaxisXS2{i}(j),plotn2{i}(j));
        end
        fprintf(fid,'};\n');
    end
end
for i = 1:nplots
    if i == 2
        fprintf(fid,'\\addplot[thick,color={rgb,255:red,%i;green,%i;blue,%i}] coordinates{\n',col0rs(i,:));
    else
        fprintf(fid,'\\addplot[color={rgb,255:red,%i;green,%i;blue,%i}] coordinates{\n',col0rs(i,:));
    end
%    fprintf(fid,'\\addplot coordinates{\n');
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
fprintf(fid,'\\legend{');
if length(varargin) == 2
    pl0ts = [plotn2;plotn];
else
    pl0ts = plotn;
end
for i = 1:length(legends)-1
    if any(pl0ts{i})
        fprintf(fid,'{%s},\n',legends{i});
    end
end
if any(pl0ts{end})
    fprintf(fid,'{%s}',legends{end});
end
fprintf(fid,'}\n',legends{end});
fprintf(fid,'\\end{semilogxaxis}\n');
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

%{
set(g,'Position',[400 200 800 350])
set(h,'XTick',[.1e-3 1 100e3 10e6],'YLim',[minval maxval],...
    'Position',[.08 .16 .89 .73]);
axis tight;

plot([1 1],[minval maxval],'--','Color',[.8 .8 .8]);
plot(100e3*[1 1],[minval maxval],'--','Color',[.8 .8 .8]);
plot([0 10e6],[1e0 1e0],'--','Color',[.8 .8 .8]);

title(mytitle);
xlabel({'Energy (eV)',xbelow});

%legend(legend1,legend2,'Location','Best')
legend(legends,'Location','Best');
%     text(.0251,xsvec(25),[sprintf('%.1f',xsvec(25))],...
%         'HorizontalAlignment','left','VerticalAlignment','bottom')

end
%}
