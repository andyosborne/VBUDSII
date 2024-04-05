classdef MuScat < NjoyParameter
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    
    end
    
    methods
        function self = MuScat(dataStruct, materialStr)
            self.nJoyType = '_3_251_';
            self.lookup = self.createLethargyAndTempLookup(dataStruct, materialStr);    
        end
    end
    
end

