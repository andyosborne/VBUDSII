function Diffusion2D1gPlots(Results,g,p)

% Results = struct('Xmesh',Xmesh,...
%                  'Ymesh',Ymesh,...
%                   'Flux',fluxOut,...
%                   'keff',keff,...
%              'powerGrid',powerEff,...
%                'runtime',runtime,...
%            'nIterations',counter-1,...
%                  'Error',ErrrVector);
%              'powerFrac',powerFrac,...
%               'gpBuffer',gpBuffer

% plotChoice 1 : flux
% plotChoice 2 : error vs iteration
% plotChoice 3 : planar flux slices

%% FINAL FLUX PLOT

if p.plotChoice(1) == 1
    Xmesh = Results.Xmesh;
    Ymesh = Results.Ymesh;
    fluxOut = Results.Flux;
    
    % Plots!
    hFlux = figure;
    % gFlux = axes;
    colormap(copper);
    
    surf(Xmesh,Ymesh,fluxOut,'EdgeColor','none');
    axis equal xy;
    view(0,90);
    
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
    
    if p.plotSave(1) == 1 
        print([p.plotName '_error.png'],'-dpng','-r300')
    end
end

if p.plotChoice(3) == 1
    Xmesh = Results.Xmesh;
    Ymesh = Results.Ymesh;
    Fluxor = Results.Flux;
    
    [sizeI sizeJ] = size(Fluxor);
    colorstr = {'r','g','b'};
    
    % width and length of reactor in physical units
    Widthby2 = 1/2*Xmesh(1,end);
    Lengthby2 = 1/2*Ymesh(end,1);
    
    Ray = [floor(.90*Widthby2); 0];
    
    figure;
    % colormap(copper);
    hold on;
    
    % Angles = pi/2*(0:3);
    Angles = pi/4*(0:7);
    Angles = pi*[0 1/8 1/4 1/2 3/4];  % <<<======== choose what slices to plot!!!
    for idx = Angles
        
        Ray2 = [cos(idx) -sin(idx); sin(idx) cos(idx)]*Ray; % rotations
        
        X = linspace(0,Ray2(1),100);
        Y = linspace(0,Ray2(2),100);
        
        R = sqrt(X.^2 + Y.^2);      % construct abscissa (ordinate?)
        
        for groupidx = 1:p.nGroups  % interpolate
            Fluxthis = interp2(Xmesh,Ymesh,Fluxor(:,:,groupidx),X+Widthby2,(Y+Lengthby2)');
            Fluxline = diag(Fluxthis);
            plot(R,Fluxline,colorstr{groupidx});
            %         fprintf('idx = %d and groupidx = %d \n',idx,groupidx);
            %         pause
        end
    end
    ylabel('Flux (cm^{-2} s^{-1})')
    xlabel('R (cm)')
%     legend('fast','resonation','thermal','Location','NorthEast');
    
    % plot inset figure showing what slices were plotted
    axes('Position',[.70 .55 .15 .15],'Layer','top','Box','off',...
        'XTick',[],'YTick',[]);
    axis equal off;
    hold on;
    plot(0,0,'ro','MarkerSize',10);
    for idx = Angles
        Ray = [cos(idx) -sin(idx); sin(idx) cos(idx)]*[0 1]';
        plot([0 Ray(1)],[0 Ray(2)],'k');
    end
    
    if p.plotSave(3) == 1 
        print([p.plotName '_slices.png'],'-dpng','-r300')
    end
end

% figure
% plot(R,thetruth-fluxOut)
