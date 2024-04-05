function out = binwise(in, varargin)

    if length(varargin) == 1
        isAbscissa = varargin{1};
    else
        isAbscissa = false;
    end


[nR nC] = size(in);
if nR == 1
    in = in';
end
if isAbscissa
    len = length(in) - 1;
    out = zeros(2*len, 1);
    out(1:2:end-1,1) = in(1:end-1);
    out(2:2:end,1) = in(2:end);
else
    len = length(in);
    out = zeros(2*len, 1);
    out(1:2:end-1,1) = in;
    out(2:2:end,1) = in;
end
end
