function Figures(StaadTable,CCC_2)

figure(1)
plot(StaadTable(:,2),StaadTable(:,3))
hold on
plot(StaadTable(:,2),StaadTable(:,4))
hold off
grid
title('Shear Force ')
xlabel('x, meter')
ylabel('FY, kN-m')

figure(2)
plot(StaadTable(:,2),StaadTable(:,5))
hold on
plot(StaadTable(:,2),StaadTable(:,6))
hold off
grid
title('Bending Moment ')
xlabel('x, meter')
ylabel('MZ, kN-m')

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
endfunction
