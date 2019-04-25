function [Sp_min,Bars_C_2,Bars_C_1,Layers_Compression,Spacing_Compression_2,Spacing_Compression_1,Bars_T_2,Bars_T_1,Layers_Tension,Spacing_Tension_2,Spacing_Tension_1]=Spacing_Tension_Compression(r,Nst_3,Nsc_5,dia_1,dia_2,B_beam,Ce,Sa)
%% spacing of bars in mm
for i=1:r  
if Nst_3(i)<=1
    SpT_provided_1=0;
  else
    SpT_provided_1=((B_beam)-(2*Ce)-(((Nst_3(i)-1)*dia_1)))/((Nst_3(i)-1));
endif

Sp_min=max([dia_1 Sa+5]);
if Sp_min<=SpT_provided_1
   SpT_Provided_1=SpT_provided_1;
   SpT_Provided_2=0;
   Layers_T=1;
   Nst_1=Nst_3(i);
   Nst_2=0;
elseif SpT_provided_1==0
   SpT_Provided_1=0;
   SpT_Provided_2=0;
   Layers_T=1;
   Nst_1=1;
   Nst_2=0;  
else
  Layers_T=2;
  Mod_T=mod(Nst_3(i),2);
    if Mod_T==0
      Nst_1=(Nst_3(i)/2);
      Nst_2=(Nst_3(i)/2);
      SpT_Provided_1=((B_beam)-(2*Ce)-((((Nst_1)-1)*dia_1)))/((Nst_1)-1);
      SpT_Provided_2=((B_beam)-(2*Ce)-((((Nst_2)-1)*dia_1)))/((Nst_2)-1);
    else
      Nst_1=ceil((Nst_3(i)/2));
      Nst_2=Nst_3(i)-ceil((Nst_3(i)/2));
      SpT_Provided_1=((B_beam)-(2*Ce)-((((Nst_1)-1)*dia_1)))/((Nst_1)-1);
      SpT_Provided_2=((B_beam)-(2*Ce)-((((Nst_2)-1)*dia_1)))/((Nst_2)-1);
    endif 
 endif
%% compression bars
if Nsc_5(i)<=1
  SpC_provided_1=0;
else
  SpC_provided_1=((B_beam)-(2*Ce)-(((Nsc_5(i)-1)*dia_2)))/((Nsc_5(i))-1);
endif

Sp_min=max([dia_1 Sa+5]);
if SpC_provided_1>=Sp_min
   SpC_Provided_1=SpC_provided_1;
   SpC_Provided_2=0;
   Layers_C=1;
   Nsc_1=Nsc_5(i);
   Nsc_2=0;
elseif SpC_provided_1==0
   SpC_Provided_1=0;
   SpC_Provided_2=0;
   Layers_C=1;
   Nsc_1=Nsc_5(i);
   Nsc_2=0;  
else
  Layers_C=2;
  Mod_C=mod(Nsc_5(i),2);
    if Mod_C==0
      Nsc_1=(Nsc_5(i)/2);
      Nsc_2=(Nsc_5(i)/2);
      SpC_Provided_1=((B_beam)-(2*Ce)-((((Nsc_1)-1)*dia_1)))/((Nsc_1)-1);
      SpC_Provided_2=((B_beam)-(2*Ce)-((((Nsc_2)-1)*dia_1)))/((Nsc_2)-1);
    else
      Nsc_1=ceil((Nsc_5(i)/2));
      Nsc_2=Nsc_5(i)-ceil((Nsc_5(i)/2));
      SpC_Provided_1=((B_beam)-(2*Ce)-((((Nsc_1)-1)*dia_1)))/((Nsc_1)-1);
      SpC_Provided_2=((B_beam)-(2*Ce)-((((Nsc_2)-1)*dia_1)))/((Nsc_2)-1);
    endif 
endif
  T(i,:)=[SpT_Provided_1];
  U(i,:)=[SpT_Provided_2];
  V(i,:)=[Layers_T];
  N1T(i,:)=[Nst_1];
  N2T(i,:)=[Nst_2];
  W(i,:)=[SpC_Provided_1];
  X(i,:)=[SpC_Provided_2];
  Y(i,:)=[Layers_C];
  N1C(i,:)=[Nsc_1];
  N2C(i,:)=[Nsc_2];
endfor

Spacing_Tension_1=[reshape(T,[],1)];
Spacing_Tension_2=[reshape(U,[],1)];
Layers_Tension=[reshape(V,[],1)];
Bars_T_1=[reshape(N1T,[],1)];
Bars_T_2=[reshape(N2T,[],1)];
Spacing_Compression_1=[reshape(W,[],1)];
Spacing_Compression_2=[reshape(X,[],1)];
Layers_Compression=[reshape(Y,[],1)];
Bars_C_1=[reshape(N1C,[],1)];
Bars_C_2=[reshape(N2C,[],1)];
Sp_min

endfunction
