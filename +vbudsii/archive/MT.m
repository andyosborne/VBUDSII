function type=MT(x)

%this program returns the MT value definitions 


        switch x
          case 1,   type='Total cross section';
          case 2,   type='* Elastic scattering cross section';
          case 3,   type='Non-elastic cross section';
          case 4,   type='(z,n) cross section';
          case 18,  type='Fission cross section';
          case 27,  type='Absorption cross section';
          case 102, type='(n,gamma) cross section';
          case 221, type='Free gas - Thermal Scattering';
          case 222, type='BOUND - Hydrogen in Water - Thermal Scattering';
          case 229, type='BOUND - Carbon in graphite - Inelastic Thermal Scattering';
          case 230, type='BOUND - Carbon in graphite - Elastic Thermal Scattering';          
          case 235, type='BOUND - Zr in ZrH - Thermal Scattering';
          case 225, type='BOUND - H in ZrH - Thermal Scattering';
          case 251, type='Mu bar - average cosine of scattering angle';   
          case 452, type='nubar - average number of neutrons released/fission';
          case '*', type='MF 6 is possible with this value - needs to be completed';
          otherwise,type='MT number is not valid or not in library.';
        end