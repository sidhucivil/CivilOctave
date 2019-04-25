clc;
clear all;
load input.mat
load data.mat
StaadTable= [csvread('Staad2.csv',1,0)];
[K]=Xd(Es,Fy)
[Mu,Vu]=Values(StaadTable);
[dia_Anchorage,Span,dia_st,B_beam,D_beam,de_beam,Ce,de_cr,C_cc]=dimensions(Span,d,Cc,Cc_Data,dia_1,dia_2,dia_3,dia_3Data,D,B,R,Shut,de_c,dia_4,dia_4Data)
[Mu_Lim]=Limiting_Moment(Fck,K,B_beam,de_beam)
[r,c] = size(StaadTable);
[Esc,Fsc]=FSC(de_cr,K,de_beam,Es,Fy,Fsc_415)
[Ast_min,Ast_max]=STEEL(B_beam,de_beam,D_beam,Fy)
fprintf('\n')
[Steel,TotalSteel_Required,TotalSteel_Provided,N_Tension,N_Compression,AreaR_Tension,AreaP_Tension,AreaP_Compression]=STEEL_R_P(dia_1,dia_2,Ast_min,Ast_max,r,Mu_Lim,Mu,Fy,Fck,de_beam,K);
[AreaP_Tension_1,Nst_3,AreaP_Tension_2]=Tension_Distribution(AreaP_Tension,dia_1,r);
[AreaP_Anchorage,AreaP_Compression_2,Nsc_5]=Compression_Distribution(AreaP_Compression,dia_2,dia_Anchorage,r);
[Sp_min,Bars_C_2,Bars_C_1,Layers_Compression,Spacing_Compression_2,Spacing_Compression_1,Bars_T_2,Bars_T_1,Layers_Tension,Spacing_Tension_2,Spacing_Tension_1]=Spacing_Tension_Compression(r,Nst_3,Nsc_5,dia_1,dia_2,B_beam,Ce,Sa);
[Tc_max]=Maximum_ShearStress(Fck,Maximum_Shear_Stress);
[Asv]=Legged_Stirrups(N_st,dia_st);
[Nominal_Shear_Stress,Design_Shear_Stress,Pt_P,Sr_R,Sv_P]=Nominal_Design_Shear(Asv,TotalSteel_Provided,B_beam,de_beam,Vu,r,Fck,Fy,Design_Shear_Strength,Design_Shear_Strength_2,Design_Shear_Strength_1);
[Tbd,Ld,CCC_1,CCC_2,CCC_3,CCC_4,CCCC,Lo]=Development_Length(Tbd_1,r,de_beam,C_cc,Tbd_switch,Fy,dia_1,Fck,Hook_Allowance,Bs,AreaP_Tension_2,B_beam,Vu);
[LD_maximum,LD_provided,Check_2]=Deflection_Check(r,AreaR_Tension,AreaP_Tension_2,B_beam,de_beam,Fy,AreaP_Compression,Modification_Factor_2,Span);
[AA_1,BB_1,R7]=Crack_Width_Check(C_cc,Mu,r,Spacing_Tension_1,Ce,dia_1,Fck,AreaP_Tension_2,B_beam,de_beam,Es,D_beam,Sigma_cbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R6=[Bars_T_1,Spacing_Tension_1,Bars_T_2,Spacing_Tension_2,Layers_Tension,Bars_C_1,Spacing_Compression_1,Bars_C_2,Spacing_Compression_2,Layers_Compression];
disp('N_T1  SP_T1  N_T2  SP_T2  Layers_T    N_CC1  SP_C1  N_C2  SP_C2  Layers_C ')
disp(num2str(R6)) 
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
disp('if Shear reinforcement not required it is given by 0 ')
disp('if Shear reinforcement required it is given by 1 ')
disp('if Shear reinforcement is not required minimum shear reinforcement is provided')
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R2=[Nominal_Shear_Stress,Pt_P,Design_Shear_Stress,Sr_R, Sv_P];
disp('    Tv      Pt_Provided       Tc     Shear reinforcement  Spacing_Shear')
disp(num2str(R2)) 
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Development length at section having maximum moment= ' num2str(Ld)])
disp('Section    MOR        Reduced_stress         LD')
disp(num2str(CCCC))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
disp('if beam stisfies the limit state of servicibility it is given by 0')
disp('if beam does not stisfies the limit state of servicibility it is given by 1')
R5=[LD_maximum,LD_provided,Check_2];
disp('Max_LD      Provided_LD        Check');
disp(num2str(R5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
disp('if beam stisfies the limit of Crack width requirement it is given by 0')
disp('if beam does not stisfies the limit of Crack width requirement it is given by 1')
disp('Section      Check');
disp(num2str(R7));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Figures(StaadTable,CCC_2)