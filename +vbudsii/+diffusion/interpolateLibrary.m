function Li = interpolateLibrary(L,p)
%INTERPOLATELIBRARY Interpolate library for a given temperature and
%Bonderenko self-shielding factor.
% Description:  The Library contains cross sections for multiple
% temperatures T and Bonderenko self-shielding factors BW. This function
% bilinearly inerpolates each cross section in the library according to a
% T and BW provided by fields in the parameter structure p.
%
% USE:  Li = interpolateLibrary(L,p)
%
% NOTES: BW values are actually exponents, but the interpolation is still
% linear between these exponent values. This is not that bad of an issue
% though, since most cross sections vary minutely with BW.
%
% EXAMPLES:
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110327 cld72   created
%   2.0      20110501 cld72   implementation after forgetting everything
%   3.0      20110519 cld72   submitting some final work this week
%
% FUTURE UPDATES:
%
% DEPENDENCIES: none
%

%     ^
%     |
%  BW | UL  UR
%     | LL  LR
%     |------->
%         T

% initialize inteprolated library
Li = L;

% grab T and BW for bilinear interpolation
T = p.T;
BW = p.BW;

% find the index of the Ts and BWs around which to interpolate
TLidx  = min( find( T - L.Ts >= 0 , 1 ,'last'), length(L.Ts)-1 );
TUidx = TLidx + 1;
BWLidx = min( find( BW - L.BWs >= 0 , 1 ,'last'), length(L.BWs)-1 );
BWUidx = BWLidx + 1;

% grab the Ts and BWs at the indices around the interpolation point
TL = L.Ts(TLidx);
TU = L.Ts(TUidx);
BWL = L.BWs(BWLidx);
BWU = L.BWs(BWUidx);

for k = L.Bins
    if k == 110 % THIS IS A HACK HERE
        for i = L.ZAIDs
            fprintf('Interpolating ZAID %d\n',i);
            for j = L.MTs

                % for brevity
                XS = L.xss(L.ZAID(i),L.MT(j)).NBins(L.Bin(k)).xs;

                % grab cross sections with which to interpolate
                if j ~= 2 % this ifstatement would be avoided with Coarse(T,BW).Vals. squeezes avoided too
                    % scattering cross section is a kernel
                    LL = squeeze(XS(TLidx,BWLidx,:));
                    LR = squeeze(XS(TUidx,BWLidx,:));
                    UL = squeeze(XS(TLidx,BWUidx,:));
                    UR = squeeze(XS(TUidx,BWUidx,:));
                else
                    % all other cross sections
                    LL = squeeze(XS(TLidx,BWLidx,:,:));
                    LR = squeeze(XS(TUidx,BWLidx,:,:));
                    UL = squeeze(XS(TLidx,BWUidx,:,:));
                    UR = squeeze(XS(TUidx,BWUidx,:,:));
                end

               % the actual bilinear interpolation
                XSout = 1/(BWU-BWL)/(TU-TL)*((BW-BWL)*((T-TL)*UR + (TU-T)*UL) ...
                    + (BWU-BW)*((T-TL)*LR + (TU-T)*LL));

                if sum(sum(sum(XSout < 0 )))
                    % it is not good if a cross seciton has
                    % a negative value, except for mubar
                    fprintf('ERROR?: Interpolated xs is negative: ZAID %d, MT %d\n',i,j);
                end

                % assign the interpolated cross section
                Li.xss(L.ZAID(i),L.MT(j)).NBins(L.Bin(k)).xs = XSout;

            end
        end
    end
end

end
