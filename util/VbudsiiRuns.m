classdef VbudsiiRuns < handle
    properties
        name = '';
        tags = {};
        fid = 0;
        params = containers.Map();
        geoms = containers.Map();
        mcnpxfluxs = containers.Map();
        Results = containers.Map();
        keffs = containers.Map();
        vbudsiifluxs = containers.Map();
        Libs = containers.Map();
        captions = containers.Map();
        savelib;
        iCool;
        iFuel;
    end
    methods
        function self = VbudsiiRuns(name, iFuel)
            self.name = name;
            self.iFuel = iFuel;
            self.iCool = 3 - self.iFuel;
        end
        function AddRun(self, tag, p, g, mcnpxflux, caption)
            self.tags{end+1} = tag;
            self.params(tag) = p;
            self.geoms(tag) = g;
            self.mcnpxfluxs(tag) = mcnpxflux;
            self.captions(tag) = caption;
        end
        function Run(self, tag, varargin)
            disp([' ************ Running ' tag]);
            [Rthis, p, g, Lthis] = vbudsii.Vbudsii(self.params(tag), ...
                                                 self.geoms(tag));
            if length(varargin) == 1
                self.savelib = varargin{1};
            else
                self.savelib = false;
            end
            if self.savelib
                save([self.name tag 'data'], 'Rthis', 'Lthis', '-v7.3');
            else
                save([self.name tag 'data'], 'Rthis', '-v7.3');
            end
        end
        function Load(self, tag)
            load([self.name tag 'data']);
            groupdef = self.params(tag).fineGroupDef;
            self.Results(tag) = Rthis;
            vbudsiiflux1 = [...
                Rthis.Region(1).Cell(self.iFuel).spectralFlux, ...
                Rthis.Region(1).Cell(self.iCool).spectralFlux, ...
                ];
            self.keffs(tag) = Rthis.Region(1).kInf;
            self.vbudsiifluxs(tag) = vbudsiiflux1 * 1e6 / ...
                sum(vbudsiiflux1(:,1) .* diff(groupdef)');
            self.mcnpxfluxs(tag) = self.mcnpxfluxs(tag) * ...
                self.GetFluxScale(tag);
            if exist('Lthis')
                self.Libs(tag) = Lthis;
            end
            % PLOT
            plotname = [self.name tag];
            vflux = self.vbudsiifluxs(tag);
            mflux = self.mcnpxfluxs(tag);
            printnuclearplot(plotname, 'flux', ...
                {['vbudsii - fuel'], ...
                    ['vbudsii - coolant'], ...
                 'mcnpx - fuel', 'mcnpx - coolant'}, ...
                groupdef, ...
                {vflux(:,1), vflux(:,2)}, ...
                'groupdefs2', groupdef, ...
                'xss2', {mflux(:,1), mflux(:,2)});
        end
        function scale = GetFluxScale(self, tag)
            groupdef = self.params(tag).fineGroupDef;
            mcnpxfluxs = self.mcnpxfluxs(tag);
            mcnpxfluxarea = sum(mcnpxfluxs(:,1) .* diff(groupdef)');
            scale = 1 / mcnpxfluxarea * 1e6;
%            scale = max(max(self.vbudsiifluxs(tag)))/ ...
%                max(max(self.mcnpxfluxs(tag)));
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
for iPlot = 1:length(self.tags)
    self.out(sprintf('%sfitzefigtwo{%s}{%s}\n\n', '\\', ...
        [self.name self.tags{iPlot}], ...
        [self.captions(self.tags{iPlot}) ...
        '. keff: ' num2str(self.keffs(self.tags{iPlot}))]));
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
