function [K]=Xd(Es,Fy)
  %% K=(Xu_max/d)
K=(0.0035*Es)/(0.0055*Es + 0.87*Fy);
endfunction
