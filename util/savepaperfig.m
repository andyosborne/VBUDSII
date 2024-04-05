function savepaperfig(fname, varargin)
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
set(0, 'DefaultAxesFontSize', 8);
set(0, 'DefaultTextFontSize', 8);
print('-depsc', fname);
%system(sprintf('ps2pdfwr -dEPSCrop %s.eps', fname));
%if isunix
%    if exist(['~/Dropbox/UTA/xspaper/figurefactory/' fname '.eps'], 'file')
%        delete(['~/Dropbox/UTA/xspaper/figurefactory/' fname '.eps']);
%    end
%    copyfile([fname '.eps'], '~/Dropbox/UTA/xspaper/figurefactory/');
%%    copyfile([fname '.pdf'], '~/Dropbox/UTA/xspaper/figurefactory');
%elseif ispc
%    copyfile([fname '.eps'], 'C:\Dropbox\UTA\xspaper\figurefactory\');
%%    copyfile([fname '.pdf'], 'C:\Dropbox\UTA\xspaper\figurefactory');
%end
end

