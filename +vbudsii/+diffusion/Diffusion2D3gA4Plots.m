function Diffusion2D3gA4Plots(Results,g,p)

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

% plotChoice 3 : top down power
% plotChoice 4 : 3d tower power

%% FINAL FLUX PLOT

if length(p.plotChoice) == 4
    
    if (p.plotChoice(3) == 1 || p.plotChoice(4) == 1)
        
        g2D = createMap4(g,p);
        
        Xmesh = Results.Xmesh;
        Ymesh = Results.Ymesh;
        Fluxor = Results.Flux;
        Assy = g2D.assy;
        PowerEff = Results.powerGrid;
        
        ReactorWidth = g2D.edgePhys;
        ReactorLength = g2D.edgePhys;
        
        % calculate normalized average power from each assy (from all 3 groups)
        % go through each assembly.
        maxPower = -inf;
        for q = 1:length(Assy)
            %PowerGrid = sum(PowerEff(Assy(q).iIdxr,Assy(q).jIdxr,:),3);
            PowerGrid = PowerEff(Assy(q).iIdxr,Assy(q).jIdxr);
            Power = mean(mean(PowerGrid)); % averaging!
            Assy(q).power = Power;
            
            if Power > maxPower
                maxPower = Power;
            end
        end
        
        if p.plotChoice(3) == 1
            % top-down plot of assy's, color indicates relative power level
            h3 = figure;
            g3 = axes;
            axis equal off;
            hold on;
            
            plot(ReactorWidth*[0 1 1 0 0],ReactorLength*[0 0 1 1 0],'k')
            
            for q = 1:length(Assy)
                fillAssy(h3,Assy(q),Assy(q).power/maxPower);
            end
            
            if p.plotSave(3) == 1
                print([p.plotName '_powerdown.png'],'-dpng','-r300')
            end
        end
        
        if p.plotChoice(4) == 1
            % 3D plot of assy power levels
            h4 = figure;
            g4 = axes;
            axis equal off;
            hold on;
            view(3);
            
            fill(ReactorWidth*[0 1 1 0 0],ReactorLength*[0 0 1 1 0],[.8 .8 .8])
            
            for q = 1:length(Assy)
                fill3Assy(h4,Assy(q),Assy(q).power/maxPower);
            end
            
            if p.plotSave(4) == 1
                print([p.plotName '_power3d.png'],'-dpng','-r300')
            end
        end
    end
end
end

function drawAssy(fhandle,Assy)

global AssySide

figure(fhandle)
plot(Assy.X+AssySide*[0 1 1 0 0],Assy.Y+AssySide*[0 0 1 1 0],'k')
    
end

function fillAssy(fhandle,Assy,color)

global AssySide

figure(fhandle)
fill(Assy.X+AssySide*[0 1 1 0],Assy.Y+AssySide*[0 0 1 1],color*[1 1 1])
    
end

function fill3Assy(fhandle,Assy,powah)

global AssySide maxPowah

m = .2;
figure(fhandle)
% X = Assy.X + AssySide*[ [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 0 0 0]', [0 1 1 0 0]' ];
% Y = Assy.Y + AssySide*[ [0 0 0 0 0]', [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 1 1 0]' ];
% Z = 10*AssySide*powah*[ [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [1 1 1 1 1]' ];
% fill3(X,Y,Z,Assy.cellid/4*[1 1 1])


X = Assy.X+m*AssySide/2 + (1-m)*AssySide*[ [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 0 0 0]', [0 1 1 0 0]' ];
Y = Assy.Y+m*AssySide/2 + (1-m)*AssySide*[ [0 0 0 0 0]', [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 1 1 0]' ];
Z = [ [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [0 0 1 1 0]', [1 1 1 1 1]' ];
fill3(X,Y,10*AssySide*powah*Z,Assy.cellid/4*[1 1 1])

m = 0;
X = Assy.X+m*AssySide/2 + (1-m)*AssySide*[ [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 0 0 0]', [0 1 1 0 0]' ];
Y = Assy.Y+m*AssySide/2 + (1-m)*AssySide*[ [0 0 0 0 0]', [0 1 1 0 0]', [1 1 1 1 1]', [1 0 0 1 1]', [0 0 1 1 0]' ];
% plot3(X,Y,10*AssySide*Z,'Color',[.9 .9 .9]);

end