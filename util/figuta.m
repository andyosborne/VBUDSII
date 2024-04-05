function figuta(hf,fname,varargin)

if isempty(varargin)
   printdir = '/home/fitze/Dropbox/UTA/r/doc/figures'; 
else
   printdir = varargin{1};
end

print(fname, '-r300', '-depsc');

system(['ps2pdfwr -dEPSCrop ' fname '.eps']);

system(['mv ' fname '.pdf ' printdir]);
system(['rm ' fname '.eps']);

end
