clc;
clear all;

load input.mat
load data.mat
StaadTable= [csvread('Staad2.csv',1,0)];


%% K=(Xu_max/d)
K=(0.0035*Es)/(0.0055*Es + 0.87*Fy)
Span

%% effective depth of beam drop plus slab thickness
if Span<=10000
 if d==0
  d_e=(Span/20);
 else
  d_e=d; 
 endif
elseif Span>10000
 if d==0
   d_e=(Span/20)*(Span/10000);
 elseif
  d_e=d;
 endif
endif


%% Clear cover
if Cc==0
  C_c=Cc_Data;
else
  C_c=Cc;
endif

if C_c<20
  disp('Minimum clear cover should be 20 mm')
  C_c=20;
endif
if Cc==0
  C_cc=max([C_c, dia_1, dia_2]);
else
  C_cc=C_c;
endif
C_cc

%% Shear reinforcement
if dia_3==0
  dia_st=dia_3Data;
else
  dia_st=dia_3;
endif
dia_st

%% Total Depth of beam
if D==0
  D_t=(d_e+C_cc+(dia_1/2)+dia_st);
else 
  D_t=(D);
  endif
%% Width of beam
if B==0
  if R>=1.5 && R<=2
     B_t=(D_t/R);
  else
    disp('D/B should be between 1.5 and 2')
    R_c=1.5
    B_t=(D_t/R_c);
  endif
 else
  B_t=B;
endif

%% Size of shuttering plates available in mm
if B==0
[m,n]=min(abs(B_t-Shut(:,1)));
B_beam= (Shut(n));
else
B_beam=B;
endif
B_beam

if D==0
[o,p]=min(abs(D_t-Shut(:,1)));
D_beam= (Shut(p));
else
D_beam=D;
endif
D_beam



%% Effective Depth of the beam
de_beam=(D_beam-C_cc-dia_st-(dia_1/2))

%% Effective cover
Ce=(D_beam-de_beam)

#figure(1)
#plot(StaadTable(:,2),StaadTable(:,3))
#hold on
#plot(StaadTable(:,2),StaadTable(:,4))
#hold off
#grid
#title('Shear Force ')
#xlabel('x, meter')
#ylabel('FY, kN-m')

figure(2)
plot(StaadTable(:,2),StaadTable(:,5))
grid
title('Bending Moment ')
xlabel('x, meter')
ylabel('MZ, kN-m')

Mu=StaadTable(:,5);
Vu=StaadTable(:,3);
Mu_Lim=(0.36*Fck*K*B_beam*de_beam*de_beam)*(1-0.416*K)/1000000

fprintf('\n');
[r,c] = size(StaadTable);

%% Effective depth of compression reinforcement
if de_c==0
  de_cr=(D_beam-de_beam);
else
  de_cr=de_c;
endif
de_cr

Esc=(0.0035)*(1-((de_cr)/(K*de_beam)))
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
    elseif
      Fsc=interp1(Fsc_500(:,1),Fsc_500(:,2),Esc)
    endif
endif
Fsc


%% Min and Max area of steel required in mm2
Ast_min=(0.85*B_beam*de_beam)/(Fy)
Ast_max=(0.04*B_beam*D_beam)
fprintf('\n')

%% Anchorage reinforcement 
if dia_4==0
  dia_Anchorage=dia_4Data;
else
  dia_Anchorage=dia_4;
endif

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

  



%%% Distribution of Tension steel
Max_Steel_Provided_Tension=max(TotalSteel_Provided)
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

R=[StaadTable];
S=[Steel];
disp('Member   Distance        Fy MAX           FY MIN          MZ MAX        MZ MIN')
disp(num2str(R))
fprintf('\n')
disp('   Mu_1            Ast_1              Mu_2          Ast_2          Asc')
disp(num2str(S))
fprintf('\n')
dia_1
dia_2
dia_Anchorage


fprintf('\n')
R1=[TotalSteel_Required,Nst_3,AreaP_Tension_2',Nsc_5',AreaP_Compression_2',AreaP_Anchorage];
disp('Steel_Required         NT     Steel_P_T              NC     Steel_P_C        Anchorage steel')
disp(num2str(R1)) 
fprintf('\n')

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
  SpC_provided_1=((B_beam)-(2*Ce)-(((Nsc_5(i)-1)*dia_1)))/((Nsc_5(i))-1);
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



R6=[Bars_T_1,Spacing_Tension_1,Bars_T_2,Spacing_Tension_2,Layers_Tension,Bars_C_1,Spacing_Compression_1,Bars_C_2,Spacing_Compression_2,Layers_Compression];
disp('N_T1  SP_T1  N_T2  SP_T2  Layers_T    N_CC1  SP_C1  N_C2  SP_C2  Layers_C ')
disp(num2str(R6)) 
fprintf('\n')








%% Tc_max(maximum shear stress in concrete) N/mm2  
 if Fck<40
    Tc_max=interp1(Maximum_Shear_Stress(:,1),Maximum_Shear_Stress(:,2),Fck);
 elseif Fck>=40
    Tc_max=4;
 endif
    Tc_max
    
%% Shear reinforcement

if N_st==2 
Asv= (N_st*pi*dia_st*dia_st)/(4)
elseif N_st==4 
Asv= (N_st*pi*dia_st*dia_st)/(4)
elseif N_st==6 
Asv= (N_st*pi*dia_st*dia_st)/(4)
else
disp("Provide 2, 4 and 6 legged only")
endif
N_st 


          ### Using X mm Y legged vertical striups


fprintf('\n')
disp('if Shear reinforcement not required it is given by 0 ')
disp('if Shear reinforcement required it is given by 1 ')
disp('if Shear reinforcement is not required minimum shear reinforcement is provided')
fprintf('\n')
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


R2=[Nominal_Shear_Stress,Pt_P,Design_Shear_Stress,Sr_R, Sv_P];
disp('    Tv      Pt_Provided       Tc     Shear reinforcement  Spacing_Shear')
disp(num2str(R2)) 
fprintf('\n')



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
disp(['Development length at section having maximum moment= ' num2str(Ld)])
disp('Section    MOR        Reduced_stress         LD')
disp(num2str(CCCC))


figure(3)
plot(StaadTable(:,2),StaadTable(:,5))
hold on
plot(StaadTable(:,2),StaadTable(:,6))
plot(StaadTable(:,2),CCC_2)
hold off
grid
title('Bending Moment ')
xlabel('x, meter')
ylabel('MZ, kN-m')


%% Deflection check
fprintf('\n')
disp('if beam stisfies the limit state of servicibility it is given by 0')
disp('if beam does not stisfies the limit state of servicibility it is given by 1')
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
R5=[LD_maximum,LD_provided,Check_2];
disp('Max_LD      Provided_LD        Check');
disp(num2str(R5));


%% Crack width check
fprintf('\n')
disp('if beam stisfies the limit of Crack width requirement it is given by 0')
disp('if beam does not stisfies the limit of Crack width requirement it is given by 1')
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
disp('Section      Check');
disp(num2str(R7));

save a.text