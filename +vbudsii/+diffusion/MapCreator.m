function regionMap = MapCreator

% creates a figure window and allows ginput of what cell goes where
% it knows how many cells/regions there are, and allows you to click where
% cell 1 is, then where cell 2 is, then where cell 3 is. then the
% reflector. it knows how fine we want the grid to be, and accordingly puts
% the right cell identifier to the correct clicked area of the reactor
% core.


h1 = figure(1);
g1 = axes;
axis equal;
hold on;

% plot reactor boundaries
side = 200;
ReactorWidth = side;
ReactorLength = side;

sidebuff = .1;


plot(ReactorWidth*[0 1 1 0 0],ReactorLength*[0 0 1 1 0],'k')
set(g1,'XLim',ReactorWidth*[-sidebuff sidebuff+1],'YLim',ReactorLength*[-sidebuff sidebuff+1])
regionMap = 1;

for j = 0:9
    x = 185-15*j;
    
    for i = 0:(9-j)
        y = 185-15*i;
        drawAssy(h1,[x y])
    end
    
end

end

function drawAssy(fhandle,lrcorner)

figure(fhandle)
AssySide = 15;
plot(lrcorner(1)+AssySide*[0 1 1 0 0],lrcorner(2)+AssySide*[0 0 1 1 0],'k')
    
end