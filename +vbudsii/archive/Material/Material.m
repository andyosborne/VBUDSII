classdef Material < handle
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name;
        crossSections;
        density;
        numberDensity;
        elScatKernel; % seems there is only a kernel for elastic
        inelScatKernel; % not used atm but maybe a nice to have in the future
        mu; % av cosine of scat angle
    end
    
    methods
        function self = Material(nameStr, dataStruct, N, rho)
            self.name = nameStr;
            self.numberDensity = N;
            self.density = rho;
            self.crossSections = self.loadCrossSections(dataStruct); %some cell struct
            self.loadScatteringKernels(dataStruct);
            self.mu = self.loadMu(dataStruct);
        end
        
        function xs = loadCrossSections(self, dataStruct)
            types = {'total' 'nGamma' 'inelastic' 'fission'};
            %types = {'elastic'};
            for i=1:length(types)
                crossSection = CrossSection(types{i}, dataStruct, self.name);
                crossSectionContainer{i} = crossSection;
            end
            xs = containers.Map(types, crossSectionContainer);
        end
        
        function loadScatteringKernels(self, dataStruct)
            self.elScatKernel = ScatteringKernel(dataStruct, self.name);
            % self.inelScatKernel = ...
        end
         
        function m = loadMu(self, dataStruct)
            m = MuScat(dataStruct, self.name);
        end
        
        function sk = getElScatKernel(self, temperature)
            sk = self.elScatKernel.lookup(temperature);
        end
       
        function m = getMuScat(self, temperature)
            m = self.mu.lookup(temperature);
        end
        
    end

    
end

