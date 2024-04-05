classdef ScatteringKernel < NjoyParameter
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function self = ScatteringKernel(dataStruct, materialStr)
            self.nJoyType = '_0_2_';
            self.lookup = self.createLethargyAndTempLookup(dataStruct, materialStr);
        end
        
        function lookupTable = createLethargyAndTempLookup(self, dataStruct, materialStr)
            for i=1:length(self.temps)
                fieldName = self.getFieldName(materialStr, self.temps{i});
                matrices{i,:,:} = self.getParamObject(dataStruct, fieldName);
            end
            lookupTable = containers.Map(self.temps, matrices);
        end
    end
    
end

