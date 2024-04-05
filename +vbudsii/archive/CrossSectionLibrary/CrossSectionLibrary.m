classdef CrossSectionLibrary < handle
    %CROSSSECTIONLIBRARY Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties
        materials;
        tempFast;
    end
    
    
    methods
        function self = CrossSectionLibrary(matFileString)
            inStruct = load(matFileString);
            self.materials = self.loadMaterials(inStruct); % some cell struct
        end
        
        function mats = loadMaterials(self, inStruct)
            names = {'U235' 'U238' 'O16' 'H2O' 'H1'};
            
            %names = {'U235'};      
            N = [0.0003867 0.007347 0.015470 0.014500 0];
            for i=1:length(names)
                matsContainer{i} = Material(names{i}, inStruct, N(i), 0);
            end
            mats = containers.Map(names, matsContainer);
        end
        
        function xss = xsSpectrum(self, materialStr, xsTypeStr, temperature)
            material = self.materials(materialStr);
            xsType = material.crossSections(xsTypeStr);
            xss = xsType.lookup(temperature);
        end
        
        function sk = elScatKernel(self, materialStr, temperature)
            material = self.materials(materialStr);
            sk = material.getElScatKernel(temperature);
        end
        
        function mu = muScat(self, materialStr, temperature)
            material = self.materials(materialStr);
            mu = material.getMuScat(temperature);
        end
    end
    
    methods(Static)
        function library = bootStrap()
            library = CrossSectionLibrary('water');
        end
        
        function library = soleInstance()
            persistent uniqueInstance;
            if isempty(uniqueInstance)
               library = CrossSectionLibrary('water');
               uniqueInstance = library;
            else
               library = uniqueInstance;
            end
        end
    end
    
end

