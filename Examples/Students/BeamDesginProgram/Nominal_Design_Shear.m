function [Nominal_Shear_Stress,Design_Shear_Stress,Pt_P,Sr_R,Sv_P]=Nominal_Design_Shear(Asv,TotalSteel_Provided,B_beam,de_beam,Vu,r,Fck,Fy,Design_Shear_Strength,Design_Shear_Strength_2,Design_Shear_Strength_1)
  %% Check for shear
for i=1:r
  %%Pt_Provided in %
  Pt_Provided=(100*TotalSteel_Provided(i))/(B_beam*de_beam);
   Tv=(Vu(i)*1000)/(B_beam*de_beam);
 
%% Design shear strength Tc in N/mm2
if Pt_Provided>=0.15
   if Fck>=40
      Tc=interp1(Design_Shear_Strength_1(:,1),Design_Shear_Strength_1(:,7),Pt_Provided);
   elseif Fck>=15
      Tc = interp2(Design_Shear_Strength,Design_Shear_Strength,Design_Shear_Strength,Fck,Pt_Provided);
   else
      Tc=disp("Value_NA");
   endif
endif

if Pt_Provided<0.15
   if Fck>=40
      Tc=0.3;
   elseif Fck>=15
      Tc=interp1(Design_Shear_Strength_2(:,1),Design_Shear_Strength_2(:,2),Fck); 
   else
      Tc=disp("Value_NA");
   endif
endif


if Pt_Provided>=3.00
    if Fck>=40
       Tc=1.01;
    elseif Fck>=15
       Tc=interp1(Design_Shear_Strength_2(:,1),Design_Shear_Strength_2(:,3),Fck); 
    else
       Tc=disp("Value_NA");
    endif
endif


if Tc>Tv
  Disp_1=0;
  Sv_1=(0.87*Fy*Asv)/(0.4*B_beam);
  Sv_2=(0.75*de_beam);
  Sv_3=(300);
  Sv=min([Sv_1 Sv_2 Sv_3]);
elseif Tv>=Tc
  Disp_1=1; 
  Tus=(Tv-Tc);
  Vus=(Tus*B_beam*de_beam);
  Sv_1P=(0.87*Fy*Asv*de_beam)/(Vus);
  Sv_1=(0.87*Fy*Asv)/(0.4*B_beam);
  Sv_2=(0.75*de_beam);
  Sv_3=(300);
  Sv=min([Sv_1P Sv_1 Sv_2 Sv_3]);
endif
E(i,:)=[Tv];
F(i,:)=[Tc];
G(i,:)=[Pt_Provided];
H(i,:)=[Disp_1];
I(i,:)=[Sv];
 endfor
Nominal_Shear_Stress=[reshape(E,[],1)];
Design_Shear_Stress=[reshape(F,[],1)];
Pt_P=[reshape(G,[],1)];
Sr_R=[reshape(H,[],1)];
Sv_P=[reshape(I,[],1)];


endfunction
