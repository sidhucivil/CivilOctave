function [AA_1,BB_1,R7]=Crack_Width_Check(C_cc,Mu,r,Spacing_Tension_1,Ce,dia_1,Fck,AreaP_Tension_2,B_beam,de_beam,Es,D_beam,Sigma_cbc)
%% Crack width check
for i=1:r
  Section=i;
Acr=sqrt(((Spacing_Tension_1(i)/2)*(Spacing_Tension_1(i)/2))+(Ce*Ce))-(dia_1/2);
Sigma_CBC= interp1(Sigma_cbc(:,1),Sigma_cbc(:,2),Fck);
Modular_ratio=280/(3*Sigma_CBC);
PtP_Provided=(100*AreaP_Tension_2'(i))/(B_beam*de_beam);
Rho=(PtP_Provided/100);
k=((2*(Rho*Modular_ratio))+(Rho*Modular_ratio*Rho*Modular_ratio))-(Rho*Modular_ratio);
x=(k*de_beam);
I_cr=((B_beam*x*x*x)/(3))+(Modular_ratio*(AreaP_Tension_2'(i))*(de_beam-x)*(de_beam-x))+((((AreaP_Tension_2'(i))*((1.5*Modular_ratio)-1))*(x-Ce)*(x-Ce)));
Fst=(((Modular_ratio*Mu(i)*1000000)*(de_beam-x))/(I_cr));
Strain_m=((D_beam-x)/(Es*(de_beam-x)))*((Fst)-(((B_beam)*(D_beam-x))/(3*AreaP_Tension_2'(i))));
W_cr=abs((3*Acr*Strain_m)/(1+(2*((Acr-C_cc)/(D_beam-x)))));
if W_cr<=0.3
  Disp_3=0;
else
  Disp_3=1;
endif
    AA(i,:)=[Section];   
    BB(i,:)=[Disp_3]; 
endfor
AA_1=[reshape(AA,[],1)];
BB_1=[reshape(BB,[],1)];
R7=[AA_1,BB_1];
endfunction
