function savetonotebookfigs(fname, varargin)

    if isempty(varargin) || strcmp(varargin{1}, 'x201u')
        copyfile(fname, '~/Dropbox/UTA/r/doc/notebook/figures')
    elseif strcmp(varargin{1}, 'x201w')
        copyfile(fname, 'C:\Dropbox\UTA\r\doc\notebook\figures')
    end
end
