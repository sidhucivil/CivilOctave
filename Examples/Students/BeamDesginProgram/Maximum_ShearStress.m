function [Tc_max]=Maximum_ShearStress(Fck,Maximum_Shear_Stress)
  
%% Tc_max(maximum shear stress in concrete) N/mm2  
 if Fck<40
    Tc_max=interp1(Maximum_Shear_Stress(:,1),Maximum_Shear_Stress(:,2),Fck);
 elseif Fck>=40
    Tc_max=4;
 endif
    Tc_max 
endfunction
