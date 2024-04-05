function out = sauerrational(x, c)
    out = 1 ./ ...
        (1 + x / c ) .^ c;
end
