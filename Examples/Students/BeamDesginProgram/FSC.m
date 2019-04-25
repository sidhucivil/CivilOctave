function [Esc,Fsc]=FSC(de_cr,K,de_beam,Es,Fy,Fsc_415)
  Esc=(0.0035)*(1-((de_cr)/(K*de_beam)));
if Fy<415
    if Esc<=0.00125
       Fsc=(0.87*Es*Esc);
    elseif Esc>0.00125
       Fsc=(0.87*Fy);
    endif
elseif Fy>=415 && Fy<500
    if Esc>=0.00380
      Fsc=360.9;
    else
      Fsc=interp1(Fsc_415(:,1),Fsc_415(:,2),Esc);
    endif
elseif Fy>=500
    if Esc>=0.00417
      Fsc=434.8;
    else
      Fsc=interp1(Fsc_500(:,1),Fsc_500(:,2),Esc);
    endif
endif
endfunction

