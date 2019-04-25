function [Asv]=Legged_Stirrups(N_st,dia_st)
%% Shear reinforcement
if N_st==2 
Asv= (N_st*pi*dia_st*dia_st)/(4);
elseif N_st==4 
Asv= (N_st*pi*dia_st*dia_st)/(4);
elseif N_st==6 
Asv= (N_st*pi*dia_st*dia_st)/(4);
else
disp("Provide 2, 4 and 6 legged only");
endif 
endfunction
