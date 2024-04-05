function Lt = InterpolateLibraryTemp(L,p)

if p.verbose
    disp('Entering InterpolateLibraryTemp');
end

Lt = L;

%% GEOFF'S CODE FROM SPECTRAL_MATRIX_ELEMENTS
%Determine how temperature selection fits into data
if p.temp<=L.Ts(1),
    temp=L.Ts(1);
    disp(['Warning: Temperature selection is below ' num2str(L.Ts(1))...
        ' K, using ' num2str(L.Ts(1)) ' K instead.']);
elseif p.temp>=L.Ts(end),
    temp=L.Ts(2);
    disp(['Warning: Temperature selection is above ' num2str(L.Ts(end))...
        ' K, using ' num2str(L.Ts(end)) ' K instead.']);
elseif sum(p.temp == L.Ts)==0 && find(abs(L.Ts-p.temp)<=5)
    temp=L.Ts( find(L.Ts > p.temp-5,1) );
    disp(['Since the reactor temperature ' num2str(temp) ' is within 5 K of '...
        num2str(temp) ' K we will use ' num2str(temp) ' K.']);
else
    temp = p.temp;
end
%%% slightly modified though.

TLoweridx = min( find( temp - L.Ts >= 0, 1, 'last'), length(L.Ts)-1);
TUpperidx = TLoweridx + 1;

TLower = L.Ts(TLoweridx);
TUpper = L.Ts(TUpperidx);

for z = L.ZAIDs % make this only the ZAIDs that matter.
    for m = L.MTs
        if m == 2
            XSL = squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,:,TLoweridx,:));
            XSU = squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,:,TUpperidx,:));
        elseif ( m == 18 || m == 452 || m == 9 ) && L.z(L.ZAID(z)).isFissile == 0
        else
            XSL = squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,TLoweridx,:));
            XSU = squeeze(L.z(L.ZAID(z)).m(L.MT(m)).xs(:,TUpperidx,:));
        end
        Lt.z(L.ZAID(z)).m(L.MT(m)).xs = 1/(TUpper-TLower)* ...
            ( (temp - TLower)*XSL + ...
              (TUpper - temp)*XSU );
    end
end
