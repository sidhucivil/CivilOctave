function [Tbd,Ld,CCC_1,CCC_2,CCC_3,CCC_4,CCCC,Lo]=Development_Length(Tbd_1,r,de_beam,C_cc,Tbd_switch,Fy,dia_1,Fck,Hook_Allowance,Bs,AreaP_Tension_2,B_beam,Vu)
%%%% Development length check
Tbd=interp1(Tbd_1(:,1),Tbd_1(:,2),Fck)
switch (Tbd_switch)
  case 1
        Ld=((0.87*Fy*dia_1)/(4* Tbd));
  case 2
        Ld=((0.87*Fy*dia_1)/(4* 1.6*Tbd));
  endswitch
Ld;
disp(['Development length at section having maximum moment= ' num2str(Ld)])
if Fy<415
  X1=3*dia_1;
elseif Fy>=415
 X1=5*dia_1;
endif
if Hook_Allowance==180
 ha=16*dia_1+4*dia_1;
elseif Hook_Allowance==90
 ha=8*dia_1+4*dia_1;
else 
disp('change the degree of hook allowance')
endif
Lo=(Bs/2)-X1-(C_cc)+ha;
for i=1
 AST_Provided=AreaP_Tension_2'(i);
 M1=(0.87*Fy*AST_Provided)*(de_beam-((Fy*AST_Provided)/(Fck*B_beam)))
 V1=Vu(i);
 Ld_1=(M1/(V1*1000));
 Z(i,:)=[Ld_1];
endfor
Ld__1=[reshape(Z,[],1)]+Lo;
disp(['Provided development_length= ' num2str(Ld__1)])
if Ld__1>=Ld
  disp('Beam satisfies the limit of Development length')
else
  disp('Beam does not satisfies the limit of Development length')
endif
fprintf('\n')
disp(['As per code Bars must extend beyond the face of support by a distance not less than ' num2str((Ld/3))])
Lif=Bs-C_cc;
disp(['Embedment Length available from inner face of support ' num2str((Lif))])

if Lif<=(Ld/3)
  disp(['There is a need to increase the embedded length'])
 Lifp=(Bs/2)+Lo;
  disp(['Embedment Length provided from inner face of support ' num2str((Lifp))])
else
 disp(['Embedment Length provided from inner face of support ' num2str((Lif))])
endif
for i=1:r
  Section_1=i;
  AST_Provided=AreaP_Tension_2'(i);
  Xu=(0.87*Fy*AST_Provided)/(0.36*Fck*B_beam);
  M2=(0.87*Fy*AST_Provided)*(de_beam-((Fy*AST_Provided)/(Fck*B_beam)));
  Zz=de_beam-((0.36/0.87)*Xu);
  sigma_stress=(M2)/(AST_Provided*Zz);
Lddd=((sigma_stress*dia_1)/(4* 1.6*Tbd));
CC(i,:)=[Section_1];
CC_1(i,:)=[M2/1000000];  
CC_2(i,:)= [sigma_stress];
CC_3(i,:)=[Lddd];
endfor
CCC_1=[reshape(CC,[],1)];
CCC_2=[reshape(CC_1,[],1)];
CCC_3=[reshape(CC_2,[],1)];
CCC_4=[reshape(CC_3,[],1)];
CCCC=[CCC_1,CCC_2,CCC_3,CCC_4];
endfunction
