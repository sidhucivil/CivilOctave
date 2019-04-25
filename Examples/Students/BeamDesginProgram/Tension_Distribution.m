function [AreaP_Tension_1,Nst_3,AreaP_Tension_2]=Tension_Distribution(AreaP_Tension,dia_1,r)
%%% Distribution of Tension steel
Max_Steel_Provided_Tension=max(AreaP_Tension)
MSP_11=0.75*Max_Steel_Provided_Tension;
MSP_12=0.5*Max_Steel_Provided_Tension;
for i=1:r
if MSP_11<=AreaP_Tension(i) && AreaP_Tension(i)<=Max_Steel_Provided_Tension
  AreaP_Tension_1(i)=Max_Steel_Provided_Tension;
  Nst_3(i,:)=[ceil((AreaP_Tension_1(i))/(0.25*pi*dia_1*dia_1))];
  AreaP_Tension_2(i)=(0.25*pi*dia_1*dia_1)*Nst_3(i);
elseif MSP_12<=AreaP_Tension(i) && AreaP_Tension(i)<MSP_11
  AreaP_Tension_1(i)=0.80*Max_Steel_Provided_Tension;
  Nst_3(i,:)=[ceil((AreaP_Tension_1(i))/(0.25*pi*dia_1*dia_1))];
  AreaP_Tension_2(i)=(0.25*pi*dia_1*dia_1)*Nst_3(i);
else
  AreaP_Tension_1(i)=0.5*Max_Steel_Provided_Tension;
  Nst_3(i,:)=[ceil((AreaP_Tension_1(i))/(0.25*pi*dia_1*dia_1))];
  AreaP_Tension_2(i)=(0.25*pi*dia_1*dia_1)*Nst_3(i);
endif
endfor

endfunction
