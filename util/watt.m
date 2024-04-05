function f = watt(E, varargin)
E = E*1e-6;
    if length(varargin) == 2
        a = varargin{1};
        b = varargin{2};
    else
        a = .2;
        b = 4;
    end

%E = 10.^(-4:.1:7)*1e-6;
%f = exp(-E/a) .* sinh(b*E).^(.5);
%C = (pi*b/4/a)^.5 * exp(b/4/a)/a;
%f = C * exp(-a*E).*sinh( (b*E).^.5);
C = sqrt(4/(pi*a^3 * b)) * exp(-a*b/4);
f = C *exp(-E/a) .* sinh((b*E).^(.5));

%semilogx(E,f)
%ylabel('probability (-)')
%xlabel('energy (MeV)')
end
