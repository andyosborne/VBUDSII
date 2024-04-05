classdef QuickReport < handle
    properties
        name = '';
        groupdef = [];
        mcnpxflux =[];
        plotnames = {};
        descriptions = {};
        fid = 0;
    end
    methods
        function self = QuickReport(name, groupdef, mcnpxflux)
            self.name = name;
            self.groupdef = groupdef;
            self.mcnpxflux = mcnpxflux;
            self.plotnames = {};
        end
        function fname = addplot(self, shortname, vbudsiiflux, varargin)
            if length(varargin) >= 2
                mcnpxflux = varargin{2};
            else
                mcnpxflux = self.mcnpxflux;
            end
            fname = [self.name shortname];
            self.plotnames{end+1} = fname;
            printnuclearplot(fname, 'flux', ...
                {['vbudsii ' shortname ' - fuel'], ...
                    ['vbudsii ' shortname ' - coolant'], ...
                 'mcnpx - fuel', 'mcnpx - coolant'}, ...
                self.groupdef, ...
                {vbudsiiflux(:,1), vbudsiiflux(:,2)}, ...
                'groupdefs2', self.groupdef, ...
                'xss2', {mcnpxflux(:,1), mcnpxflux(:,2)});
            if length(varargin) >= 1
                self.descriptions{end+1} = varargin{1};
            else
                self.descriptions{end+1} = '';
            end
        end
        function printreport(self)
fname = [self.name '.tex'];
self.fid = fopen(fname, 'w');
self.out('\\documentclass[letterpaper,12pt]{article}\n');
self.out('\\usepackage{graphicx}\n');
self.out('\\usepackage[left=.5in,right=.5in,top=1in,bottom=1in]{geometry}\n');
self.out('\\usepackage{caption}\n');
self.out('\\newcommand{\\fitzefigtwo}[2]{\\begin{center}\n');
self.out('\\includegraphics{#1}\\captionof{figure}{#2}\\label{fig:#1}\n');
self.out('\\end{center}}\n');
self.out('\\begin{document}\n');
self.out(sprintf('%stitle{%s}\n', '\\', self.name));
self.out('\\date{\\today}\n');
self.out('\\author{Chris Dembia}\n');
self.out('\\maketitle\n');
for iPlot = 1:length(self.plotnames)
    self.out(sprintf('%sfitzefigtwo{%s}{%s}\n\n', '\\', ...
        self.plotnames{iPlot},...
        self.descriptions{iPlot}));
end
self.out('\\end{document}\n');
fclose(self.fid);
system(sprintf('rubber --pdf %s', fname));
        end
        function out(self, string)
            fprintf(self.fid, string);
        end

    end
end
