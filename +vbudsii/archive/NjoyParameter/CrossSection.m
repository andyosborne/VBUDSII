classdef CrossSection < NjoyParameter
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type;
        typeMap = containers.Map({'total' 'nGamma' ...
            'inelastic' 'fission'},...
            {'_3_1_' '_3_102_' '_3_4_' '_3_18_'});
    end
    
    methods
        function self = CrossSection(typeStr, dataStruct, materialStr)
            self.type = typeStr;
            self.nJoyType = self.typeMap(typeStr);
            %self.nJoyType = getfield(self.typeMap, typeStr);
            self.lookup = self.createLethargyAndTempLookup(dataStruct, materialStr);
        end
        

    end
    
end

