function Lts = InterpolateLibrarySShield(L,p)

if p.verbose
    disp('Entering InterpolateLibrarySShield');
end

Lts = L;

s0 = 10; % when comparing to Geoff's code use this.

if isfield(p,'constS0') == 1
    s0 = p.constS0;
end

for z = L.ZAIDs
    for m = L.MTs
        if m == 2
            Lts.z(L.ZAID(z)).m(L.MT(m)).xs = ...
                L.z(L.ZAID(z)).m(L.MT(m)).xs(:,:,L.S0(s0));
        elseif ( m == 18 || m == 452 || m == 9 ) && ...
            L.z(L.ZAID(z)).isFissile == 0
        else
            Lts.z(L.ZAID(z)).m(L.MT(m)).xs = ...
                L.z(L.ZAID(z)).m(L.MT(m)).xs(:,L.S0(s0));
        end
    end
end
