function paperfigureprops(varargin)
if length(varargin) == 1
    doublewidth = varargin{1};
else
    doublewidth = false;
end
imageWidth = 8.5;
imageHeight = 8.0;
if doublewidth
    imageWidth = 2 * imageWidth;
end
set(gcf, 'PaperUnits', 'centimeters', ...
         'PaperPosition', [0 0 imageWidth imageHeight]);
set(0, 'DefaultAxesFontSize', 10);
set(0, 'DefaultTextFontSize', 10);

end
