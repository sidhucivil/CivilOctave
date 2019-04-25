function [Mu_Lim]=Limiting_Moment(Fck,K,B_beam,de_beam)
  Mu_Lim=(0.36*Fck*K*B_beam*de_beam*de_beam)*(1-0.416*K)/1000000
endfunction
