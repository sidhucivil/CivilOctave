function [Steel,TotalSteel_Required,TotalSteel_Provided,N_Tension,N_Compression,AreaR_Tension,AreaP_Tension,AreaP_Compression]=STEEL_R_P(dia_1,dia_2,Ast_min,Ast_max,r,Mu_Lim,Mu,Fy,Fck,de_beam,K)
for i=1:r
%% Singly
 if Mu_Lim>=Mu(i)
   Mu_1=Mu(i);
   Ast_1=(Mu(i)*1000000)/(0.87*Fy*de_beam*(1-0.416*K));
   Mu_2=0;
   Ast_2=0;
   Asc=0;
%% Doubly
 elseif Mu_Lim<Mu(i)
   Mu_1=Mu_Lim;
   Ast_1=(Mu_1*1000000)/(0.87*Fy*de_beam*(1-0.416*K));
   Mu_2=(Mu(i)-Mu_Lim);
   Ast_2=(Mu_2*1000000)/(0.87*Fy*(de_beam-de_cr));
   Asc=(Mu_2*1000000)/((Fsc-0.446*Fck)*(de_beam-de_cr));
  endif
AstR_T=Ast_1+Ast_2;
As_TC=Ast_1+Ast_2+Asc ;
if As_TC>=Ast_min && As_TC<=Ast_max
    AstR_1=As_TC;
  elseif As_TC<Ast_min
   AstR_1=Ast_min;
  elseif As_TC>Ast_max
   AstR_1=Ast_max;
  endif

A(i,:)=[Mu_1 Ast_1 Mu_2 Ast_2 Asc];
B(i,:)=[AstR_1];
J(i,:)=[AstR_T];
Nst(i,:)=[ceil((Ast_1+Ast_2)/(0.25*pi*dia_1*dia_1))];
Nsc(i,:)=[ceil((Asc/(0.25*pi*dia_2*dia_2)))];
AsTP_1=(Nst(i,:)*0.25*pi*dia_1*dia_1);
AsCP_2=(Nsc(i,:)*0.25*pi*dia_2*dia_2);
C(i,:)=[AsTP_1+AsCP_2];
C_1(i,:)=[AsTP_1];
C_2(i,:)=[AsCP_2];

endfor
Steel=[reshape(A,[],5)];
TotalSteel_Required=[reshape(B,[],1)];
TotalSteel_Provided=[reshape(C,[],1)];
N_Tension=[reshape(Nst,[],1)];
N_Compression=[reshape(Nsc,[],1)];
AreaR_Tension=[reshape(J,[],1)];
AreaP_Tension=[reshape(C_1,[],1)];
AreaP_Compression=[reshape(C_2,[],1)];