classdef NjoyParameter < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lookup;
        nJoyType;
        temps = {300 600 900 1200 1500};
    end
    
    methods
        function fName = getFieldName(self, materialStr, temp)
            fName = [materialStr self.nJoyType int2str(temp)];
        end
        
        function parameterObject = getParamObject(self, dataStruct, fieldName)
            exists = isfield(dataStruct, fieldName);
            if exists == 1
                parameterObject = getfield(dataStruct, fieldName);
            else
                parameterObject = 0;
            end
        end
        
        function lookupTable = createLethargyAndTempLookup(self, dataStruct, materialStr)
            for i=1:length(self.temps)
                fieldName = self.getFieldName(materialStr, self.temps{i});
                vectors{i,:} = self.getParamObject(dataStruct, fieldName);  
            end
            lookupTable = containers.Map(self.temps, vectors);
        end        
    end
end

