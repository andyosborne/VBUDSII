function uc = initunitcell()

uc = struct('pinPitch', [], ... % 1 x nAntipins
            'pinDiam', [], ... % 1 x nPins
            'f', [], ... % alpha, 1 x nAntipins
            'g', [], ... % beta, 1 x nPins
            'sauerConst', struct('mod', [], ... % scalar
                                 'fuel', []), ... % scalar
            'PinCellNames', {}, ... % 1 x nPins
            'AntipinCellNames', {}); % 1 x nAntipins
                                  
            
end
