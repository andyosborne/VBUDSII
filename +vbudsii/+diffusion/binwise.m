function out = binwise(in)

if (length(in) ~= 110)
    error('Only works for things that are 110 long.')
end

out(1:2:219,1) = in;
out(2:2:220,1) = in;

end