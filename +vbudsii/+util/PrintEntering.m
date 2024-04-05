function PrintEntering(p, string)

distance = 30;

if p.verbose
    disp(['  ----- ' string ': ' blanks(distance - length(string)) ...
        'entering -----'])
end

end
