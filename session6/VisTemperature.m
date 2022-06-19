function VisTemperature(formfunction, X, Y, T, l, str)



figure(2)
set(gca, 'FontSize', 18)
set(gcf,'position', [34   164   560   420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
hold off
surf(X, Y, T);

hold on
campos([58.6464  -61.7142  130.2245])
surf(X, -Y, T)
zlabel('T')
xlabel('x')
ylabel('y')
zlim([ min([min(T), 90])  max([100, max(T)])   ])
if nargin == 6
title(str)
end




 set(gcf, 'color', 'white')
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])



figure(3)
set(gca, 'FontSize', 18)
set(gcf,'position', [613   163   560   420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
 set(gcf, 'color', 'white')
hold off
plot([0 linspace(0, l, 50) l], [0 formfunction(linspace(0, 1, 50)) 0])
hold on
plot([0 linspace(0, l, 50) l], [0 -formfunction(linspace(0, 1, 50)) 0])
contour(X,Y,T)
hold on
contour(X, -Y, T)

% contour(X,-Y,T)

ylim([-5, 5])
axis equal

xlabel('x')
ylabel('y')
if nargin == 6
title(str)
end
colorbar

figure(4)
set(gca, 'FontSize', 18)
set(gcf,'position', [ 1190         162         560         420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
 set(gcf, 'color', 'white')
hold off
plot([0 linspace(0, l, 50) l], [0 formfunction(linspace(0, 1, 50)) 0])
hold on
plot([0 linspace(0, l, 50) l], [0 -formfunction(linspace(0, 1, 50)) 0])
pcolor(X,Y,T)
shading interp
pcolor(X,-Y,T)
shading interp

ylim([-5, 5])
axis equal

xlabel('x')
ylabel('y')

if nargin == 6
title(str)
end
colorbar

