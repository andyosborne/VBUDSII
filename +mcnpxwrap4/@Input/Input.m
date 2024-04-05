classdef Input < handle
    properties
        param;
        geom;
        inelastic = false;
    end
    methods
        function self = Input(vbudsiiP, vbudsiiG, inelastic)
            self.param = vbudsiiP;
            self.geom = vbudsiiG;
            self.inelastic = inelastic;
        end
        Write(self, mcnpxexecpath, ifname, runhard);
    end
end
