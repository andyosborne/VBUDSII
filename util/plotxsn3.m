function g = plotxsn3(groupdefs, xs, varargin)

myInputParser = inputParser;

defaultBinit = true;
validBinit = {true, false};

defaultUnits = 'micro';
validUnits = {'micro', 'macro'};
checkUnits = @(x) any(validatestring(x, validUnits));

defaultTitle = '';

defaultXlabelBelow = '';

p.addOptional('binit', defaultBinit)
p.addOptional('units', defaultUnits, checkUnits);
p.addOptional('title', defaultTitle, @ischar);
p.addOptional('xlabelbelow', defaultXlabelBelow, @ischar);

% legends, and multiple cross sections.
% groupdefs can be 1D or cell array.

binit,macro,mytitle,legends,xbelow,groupdefs,xss)

    nplots = length(xss);
    lengthX = zeros(nplots,1);
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
% prepare for bin-wise plotting.
plotn = cell(nplots,1);
minval = inf;
maxval = -inf;
for i = 1:nplots
    if binit == 1
        plotn{i} = zeros(lengthX(i)*2,1);
        plotn{i}(1:2:2*lengthX(i)-1) = xss{i};
        plotn{i}(2:2:2*lengthX(i)) = xss{i};
    else
        plotn{i} = xss{i};
    end

    %if max(plotn(:,i)) > maxval
    %    maxval = max(plotn{i});
    %end
    %if min(plotn(:,i)) < minval
    %    minval = min(plotn{i});
    %end
end

% plot vector cross section.
g = figure;
h = axes;
%loglog(XaxisXS*ones(1,nargin-4),plotn)
%handles = loglog(XaxisXS*ones(1,nargin-6),plotn);
%set(handles(1),'LineWidth',2);
styles = {'b','g','r','c','m','k'};
handles = zeros(nplots,1);
for i = 1:nplots
    handles(i) = loglog(XaxisXS{i},plotn{i},styles{i});
    if i == 1
        hold on;
    end
end

set(g,'Position',[400 200 800 350])
set(h,'XTick',[.1e-3 1 100e3 10e6],... %'YLim',[minval maxval],...
    'Position',[.08 .16 .89 .73]);
%axis tight;

plot([1 1],[minval maxval],'--','Color',[.8 .8 .8]);
plot(100e3*[1 1],[minval maxval],'--','Color',[.8 .8 .8]);
plot([0 10e6],[1e0 1e0],'--','Color',[.8 .8 .8]);

title(mytitle);
xlabel({'Energy (eV)',xbelow});

if macro == 1
    ylabel('\Sigma (cm^{-1})');
elseif macro == 0
    ylabel('$\sigma$ (b)','Interpreter','latex');
elseif macro == 2
    ylabel('$\phi (n cm^{-2} s^{-1})$','Interpreter','latex')
end
%legend(legend1,legend2,'Location','Best')
legend(legends,'Location','Best');
%     text(.0251,xsvec(25),[sprintf('%.1f',xsvec(25))],...
%         'HorizontalAlignment','left','VerticalAlignment','bottom')

end


% 
% plotxs(0,H1_3_102_600_10,'H1 microscopic absorption')
% close all
% plotxs(0,U235_3_18_600_10+U235_3_102_600_10,'U235 microscopic absorption')
% plotxs(0,U238_3_18_600_10+U238_3_102_600_10,'U235 microscopic absorption')
% plotxs(0,U238_3_18_600_10+U238_3_102_600_10,'U238 microscopic absorption')
% plotxs(0,O16_3_102_600_10,'U238 microscopic absorption')
% plotxs(0,O16_3_102_600_10,'O16 microscopic absorption')
% plotxs(0,sum(H1_0_2_600),'H1 microscopic scattering')
% plotxs(0,[sum(H1_0_2_600)' sum(H1_0_2_600(1:40,:))'],'H1 microscopic scattering')
% plotxs(0,sum(H1_0_2_600(1:40,:)),'H1 microscopic scattering')
% plotxs(0,sum(H1_0_2_600),'H1 microscopic scattering')
% plotxs(0,sum(H1_0_2_600),'H1 microscopic scattering, just to thermal')
% plotxs(0,sum(H1_0_2_600(1:40,:)),'H1 microscopic scattering, just to thermal')

%  
% plotxs(0,sum(H1_0_2_600(1:40,:)),'H1 microscopic scattering')
% plotxs(0,sum(H1_0_2_600),'H1 microscopic scattering')
% plotxs(0,sum(H1_0_2_600),'H1 microscopic scattering, just to thermal')
% plotxs(0,sum(H1_0_2_600(1:40,:)),'H1 microscopic scattering, just to thermal')
% plotxs(0,sum(O16_0_2_600(1:40,:)),'O16 microscopic scattering, just to thermal')
% plotxs(0,sum(O16_0_2_600),'O16 microscopic scattering')
% plotxs(0,sum(O16_0_2_600(41:90,:)),'O16 microscopic scattering')
% plotxs(0,sum(H2O_0_2_600),'H2O microscopic scattering')
% plotxs(0,sum(H2O_0_2_600(1:40,:)),'H2O microscopic scattering, just to thermal')
