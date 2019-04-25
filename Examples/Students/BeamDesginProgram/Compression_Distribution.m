function [AreaP_Anchorage,AreaP_Compression_2,Nsc_5]=Compression_Distribution(AreaP_Compression,dia_2,dia_Anchorage,r)
  %%% Distribution of Compression steel
Max_Steel_Provided_Comp=max(AreaP_Compression)
MSPC_11=0.75*Max_Steel_Provided_Comp;
MSPC_12=0.6*Max_Steel_Provided_Comp;
for i=1:r
if Max_Steel_Provided_Comp==0
   AreaP_Compression_2(i)=0;
   Nsc_5(i)=0;
  elseif MSPC_11<=AreaP_Compression(i) && AreaP_Compression(i)<=Max_Steel_Provided_Comp
  AreaP_Compression_1(i)=Max_Steel_Provided_Comp;
  Nsc_4(i,:)=[ceil((AreaP_Compression_1(i))/(0.25*pi*dia_2*dia_2))];
   if Nsc_4(i)>0 && Nsc_4(i)<=2
      Nsc_5(i)=2;
   elseif Nsc_4(i)>2
      Nsc_5(i)=Nsc_4(i);
   endif
  AreaP_Compression_2(i)=(0.25*pi*dia_2*dia_2)*Nsc_5(i);
 elseif MSPC_12<=AreaP_Compression(i) && AreaP_Compression(i)<MSPC_11
  AreaP_Compression_1(i)=0.80*Max_Steel_Provided_Comp;
  Nsc_4(i,:)=[ceil((AreaP_Compression_1(i))/(0.25*pi*dia_2*dia_2))];
  if Nsc_4(i)>0 && Nsc_4(i)<=2
      Nsc_5(i)=2;
   elseif Nsc_4(i)>2
      Nsc_5(i)=Nsc_4(i);
   endif
  AreaP_Compression_2(i)=(0.25*pi*dia_2*dia_2)*Nsc_5(i);  
 elseif AreaP_Compression(i)<MSPC_12
  AreaP_Compression_1(i)=0.6*Max_Steel_Provided_Comp;
  Nsc_4(i,:)=[ceil((AreaP_Compression_1(i))/(0.25*pi*dia_2*dia_2))];
  if Nsc_4(i)>0 && Nsc_4(i)<=2
      Nsc_5(i)=2;
   elseif Nsc_4(i)>2
      Nsc_5(i)=Nsc_4(i);
   endif
  AreaP_Compression_2(i)=(0.25*pi*dia_2*dia_2)*Nsc_5(i);
 endif 
 if AreaP_Compression_2==0
  %%provide 2 anchorage bars on the top to hold the stirrups
    As_Anchorage=2*(0.25*pi*dia_Anchorage*dia_Anchorage);
  else
    As_Anchorage=0;
 endif
   L(i,:)=[ As_Anchorage];
endfor
AreaP_Anchorage=[reshape(L,[],1)];
endfunction
