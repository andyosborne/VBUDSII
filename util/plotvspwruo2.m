function plotvspwruo2(vbudsiiflux)
load('~/Dropbox/UTA/r/code/util/mcnpxpwruo2flux'); %, 'mcnpxflux');
mcnpxflux = mcnpxflux * max(max(vbudsiiflux)) / max(max(mcnpxflux));
plot([vbudsiiflux mcnpxflux]);
legend('vbudsii fuel', 'vbudsii coolant', 'mcnpx fuel', 'mcnpx coolant');
end
