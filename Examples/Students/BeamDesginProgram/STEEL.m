function [Ast_min,Ast_max]=STEEL(B_beam,de_beam,D_beam,Fy)
  %% Min and Max area of steel required in mm2
Ast_min=(0.85*B_beam*de_beam)/(Fy);
Ast_max=(0.04*B_beam*D_beam);
endfunction
