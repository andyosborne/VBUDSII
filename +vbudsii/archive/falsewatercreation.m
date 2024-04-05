function [output]=watercreation
load data


MT=[0 3 6];
MF=[1 2 102 221 222 251];
TEMP=[300 600 900 1200 1500]
TEMPw=[296 600 1000]
MFMT={'_0_2_', '_3_1_', '_3_2_', '_6_2_', '_3_102_', '_3_221_', '_6_221_', '_3_222_', '_6_222_'}; % this is a cell array, use char(MFMT(2)) to make a string.
%reorder for better operation (collect in like operations


%here we will just create water at 600 degrees

%for n=1:1:

% eval([  'H2O' char(MFMT(1)) =  
     
 
 
 
 
 %this whole thing would be easier if I restricted by 221 and 222 to the
 %thermal regions.  Perhaps doing thermal separately is a good idea.
 
 
 
 
 
 %2 222 221 Scatter cross section
 
 HinH2O= H1_0_2_600;
 HinH2O(1:4,1:4)= H1_6_221_600(1:4,1:4);
 
 H2O_0_2_600=2*HinH2O+O16_0_2_600;
 
 %102 Radiative capture cross section
 
 HinH2O= H1_3_102_600;
 HinH2O(1:4)= H1_3_102_600(1:4);
 
 H2O_3_102_600=2*HinH2O+O16_3_102_600;
 
 %251 % same as above but 251 in place of 102
 
 HinH2O= H1_3_251_600;
 HinH2O(1:4)= H1_3_251_600(1:4);
 
 H2O_3_251_600=2*HinH2O+O16_3_251_600;
 
 
 %1 % same as above but 1 in place of 102
 
 HinH2O= H1_3_1_600;
 HinH2O(1:4)= H1_3_1_600(1:4);
 
 H2O_3_1_600=2*HinH2O+O16_3_1_600;
 
  H2O_3_4_600=O16_3_4_600;
 
 save falsewater
 %falsewater is water without the chemical bonding energy being used.  