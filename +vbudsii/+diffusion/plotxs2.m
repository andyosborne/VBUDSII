function g = plotxs2(macro,xsvec1,xsvec2,mytitle,legend1,legend2,xbelow)

load FLUXDATA
XaxisXS = Xaxis(end:-1:1);
% prepare for bin-wise plotting.
plot1 = binwise(xsvec1);
plot2 = binwise(xsvec2);

    % plot vector cross section.
    g = figure;
    h = axes;
    loglog(XaxisXS,plot1,'r')
    hold on;
    loglog(XaxisXS,plot2,'b-')

    set(g,'Position',[400 200 800 350])
    set(h,'XTick',[.1e-3 1 100e3 10e6],'YLim',[min(min(plot1,plot2)) max(max(plot1,plot2))],'Position',[.08 .16 .89 .73])
    axis tight
    
    plot([1 1],[min(min(plot1,plot2)) max(max(plot1,plot2))],'--','Color',[.8 .8 .8])
    plot(100e3*[1 1],[min(min(plot1,plot2)) max(max(plot1,plot2))],'--','Color',[.8 .8 .8])
    
    title(mytitle)
    xlabel({'Energy (eV)',xbelow})

    if macro == 1
        ylabel('\Sigma (cm^{-1})')
    elseif macro == 0
        ylabel('\sigma (b)')
    end
    legend(legend1,legend2,'Location','Best')
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
