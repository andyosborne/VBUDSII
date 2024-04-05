function out = carlvikrational(x)

    b1 = 2;
    b2 = -1;
    a1 = 2;
    a2 = 3;
    out = 1 - x .* ( b1 ./ (x + a1) + b2 ./ (x + a2) );

end
