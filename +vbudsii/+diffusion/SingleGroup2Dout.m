function [fluxIntOut] = SingleGroup2Dout(source,fluxIntIn,A,BoundCond) 
%SINGLEGROUP2DOUT Solves one-group diffusion on a 2D uniform grid, given a
%source, grid coefficients, and boundary conditions
% Description:  Provides single-group flux based on a neutron source term,
% a flux guess, grid point coefficients, and boundary conditions. Flux
% output is internal flux, which is the flux except at the outer-most
% square of grid points. Those outer grid points are used to specify
% boundary conditions.
%
% USE:  [fluxIntOut] = SingleGroup2Dout(source,fluxIntIn,A,BoundCond) 
%
% NOTES: Implementation of Stacey 2001 (p83). Diffusion2D is the only
% function that calls this function.
%
% EXAMPLES: 
%
% MAJOR UPDATES:
%   version  date     NetID   description
%   1.0      20110520 cld72   cleaned up and formatted
%              
% FUTURE UPDATES:
%
% DEPENDENCIES: none
%

global numNodeI numNodeJ

% initialize fluxes

fluxA = zeros(numNodeI,numNodeJ);
fluxA(2:numNodeI-1,2:numNodeJ-1) = fluxIntIn; % flux from "last iteration", or guess
fluxB = fluxA; % new flux
%% Process Boundary conditions

% ZERO FLUX

% left
if BoundCond(1) == 0
    fluxA(1:numNodeI,1) = zeros(numNodeI,1);
end

% bottom (physically)
if BoundCond(2) == 0
    fluxA(1,1:numNodeJ) = zeros(1,numNodeJ);
end

% right
if BoundCond(3) == 0
    fluxA(1:numNodeI,numNodeJ) = zeros(numNodeI,1);
end

% top (physically)
if BoundCond(4) == 0
    fluxA(numNodeI,1:numNodeJ) = zeros(1,numNodeJ);
end

% symmetry
if BoundCond(1) == 1
    fluxA(2:numNodeI-1,1) = fluxA(2:numNodeI-1,2);
%     fluxA(1,1) = fluxA(2,2);
%     fluxA(numNodeI,1) = fluxA(numNodeI-1,2);
end

if BoundCond(2) == 1
    fluxA(1,2:numNodeJ-1) = fluxA(2,2:numNodeJ-1);
%     fluxA(1,1) = fluxA(2,2);
%     fluxA(1,numNodeJ) = fluxA(2,numNodeJ-1);    
end

if BoundCond(3) == 1
    fluxA(2:numNodeI-1,numNodeJ) = fluxA(2:numNodeI-1,numNodeJ-1);
%     fluxA(1,1) = fluxA(2,2);
%     fluxA(numNodeI,numNodeJ) = fluxA(numNodeI-1,numNodeJ-1);
end

if BoundCond(4) == 1
    fluxA(numNodeI,2:numNodeJ-1) = fluxA(numNodeI-1,2:numNodeJ-1);
%     fluxA(numNodeI,1) = fluxA(numNodeI-1,2);
%     fluxA(numNodeI,numNodeJ) = fluxA(numNodeI-1,numNodeJ-1);
end

% THE SYMMETRY CONDITION MUST BE CONTINUALLY UPDATED INSIDE THE LOOP

% set up Gauss-Seidel relaxation
maxIterations = 2000;
isErrorful = 1;
itcounter = 1;
thresholdp = 1e-2;
ErrrVector = [];
while itcounter <= maxIterations && isErrorful

    p = 1;
    for j = 2:numNodeJ-1
        for i = 2:numNodeI-1

            % this calculation requires knowledge of boundary fluxes.
            I = i-1;
            J = j-1;
            fluxB(i,j) = source(I,J)...
                       - A(I,J,2)*fluxB(i-1,j) - A(I,J,4)*fluxB(i,j-1)...
                       - A(I,J,3)*fluxA(i+1,j) - A(I,J,5)*fluxA(i,j+1);
                          
            fluxB(i,j) = fluxB(i,j)/A(I,J,1);
  
            p = p + 1;
            
            % change A according to boundary conditions
            
        end
    end
    
    % take care of boundary conditions; assign to boundary the correct flux.
   
    % control error from old flux and just-calculated flux
    isErrorful = norm(fluxA-fluxB)/norm(fluxA) > thresholdp;
    
    errr = norm(fluxA-fluxB);
    ErrrVector = [ErrrVector errr];
    
    itcounter = itcounter + 1;
    
    % update boundary if boundcond == 1    
    % SYMMETRY
    if BoundCond(1) == 1
        fluxB(2:numNodeI-1,1) = fluxB(2:numNodeI-1,2);
    end
    
    if BoundCond(2) == 1
        fluxB(1,2:numNodeJ-1) = fluxB(2,2:numNodeJ-1);
    end
    
    if BoundCond(3) == 1
        fluxB(2:numNodeI-1,numNodeJ) = fluxB(2:numNodeI-1,numNodeJ-1);
    end
    
    if BoundCond(4) == 1
        fluxB(numNodeI,2:numNodeJ-1) = fluxB(numNodeI-1,2:numNodeJ-1);
    end
    
    % current flux becomes old source       
    fluxA = fluxB;
    
end

% if convergence does not occur in a reasonable number of iterations
if itcounter > maxIterations
    error('SingleGroup2D has not converged. put more details in this error msg')
end

% only send out internal flux
fluxIntOut = fluxB(2:numNodeI-1,2:numNodeJ-1);

end