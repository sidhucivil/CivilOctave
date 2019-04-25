function [dia_Anchorage,Span,dia_st,B_beam,D_beam,de_beam,Ce,de_cr,C_cc]=dimensions(Span,d,Cc,Cc_Data,dia_1,dia_2,dia_3,dia_3Data,D,B,R,Shut,de_c,dia_4,dia_4Data)
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
 else
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
C_cc;
%% Shear reinforcement
if dia_3==0
  dia_st=dia_3Data;
else
  dia_st=dia_3;
endif
dia_st;

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
    R_c=1.5;
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
B_beam;
if D==0
[o,p]=min(abs(D_t-Shut(:,1)));
D_beam= (Shut(p));
else
D_beam=D;
endif
D_beam;

%% Effective Depth of the beam
de_beam=(D_beam-C_cc-dia_st-(dia_1/2));
%% Effective cover
Ce=(D_beam-de_beam);
%% Effective depth of compression reinforcement
if de_c==0
  de_cr=(D_beam-de_beam);
else
  de_cr=de_c;
endif
de_cr;
%% Anchorage reinforcement 
if dia_4==0
  dia_Anchorage=dia_4Data;
else
  dia_Anchorage=dia_4;
endif
endfunction
