function Diffusion1D3gPlots(Results,g,p)

%     Results = struct('R',R,...
%                   'Flux',fluxOut,...
%                   'keff',keff,...
%              'powerFrac',powerFrac,...
%                'runtime',runtime,...
%            'nIterations',counter-1,...
%                  'Error',ErrrVector);

% plotChoice 1 : flux
% plotChoice 2 : error vs iteration
%% FINAL FLUX PLOT
if p.plotChoice(1) == 1
    R = Results.R;
    fluxOut = Results.Flux;
    
    hfinalflux = figure;
    gfinalflux = axes;
    hold on;
    plot(R*[1 1 1],fluxOut)
    legend('fast','resonance','thermal') % EDIT resonation
    xlabel('R (cm)');
    % xlabel({'R (cm)','',['Number of iterations: ' num2str(counter-1)],['Runtime: ' num2str(runtime) ' seconds']})
    ylabel('\phi (# cm^{-2} s^{-1})')
    
    bottoM = min(min(fluxOut));
    toP = max(max(fluxOut));
    % for i = 1:length(R)
    %     plot(R(i)*[1 1],[min(bottoM,0) 1.1*toP],':','Color',[.8 .8 .8])
    % end
    
    % for i = 2:4
    %     plot(campR(i)*[1 1],[min(bottoM,0) 1.1*toP],'k')%'LineWidth',1.5)
    % end
    
    set(gfinalflux,'XLim',[0 R(end)],'YLim',[0 1.1*toP],'Box','on')
    
    if p.plotSave(1) == 1
        print([p.plotName '_flux.png'],'-dpng','-r300')
    end
end

if p.plotChoice(2) == 1
    ErrrVector = Results.Error;
    hError = figure;
    gError = axes;
    semilogy(ErrrVector)
    xlabel('iteration (-)')
    ylabel('approximate error (-)')
    
    if p.plotSave(2) == 1
        print([p.plotName '_error.png'],'-dpng','-r300')
    end
    
end
