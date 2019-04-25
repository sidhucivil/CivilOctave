function [LD_maximum,LD_provided,Check_2]=Deflection_Check(r,AreaR_Tension,AreaP_Tension_2,B_beam,de_beam,Fy,AreaP_Compression,Modification_Factor_2,Span)
%% Deflection check
for i=1:r
   PtP_Provided=(100*AreaP_Tension_2'(i))/(B_beam*de_beam);
   Fs=(0.58*Fy)*(AreaR_Tension(i)/AreaP_Tension_2'(i));
   Kt=(1/((0.225)+(0.00322*Fs)-(0.625*log10(1/PtP_Provided))));
   PtP_Compression=(AreaP_Compression(i)/(B_beam*de_beam))*100;
   Kc=interp1(Modification_Factor_2(:,1),Modification_Factor_2(:,2),PtP_Compression);
   Kf=1;  
   %% L/Dmax=LD
   LD_max=20*Kt*Kc*Kf;
   LD_Provided=(Span)/(de_beam);
   if LD_max>=LD_Provided
     Disp_2=0;
   elseif LD_max<LD_Provided
     Disp_2=1;
   endif
    O(i,:)=[LD_max];   
    P(i,:)=[LD_Provided]; 
    Q(i,:)=[Disp_2];
endfor
LD_maximum=[reshape(O,[],1)];
LD_provided=[reshape(P,[],1)];
Check_2=[reshape(Q,[],1)];
endfunction
