% function VisTemperature_transient(formfunction, X, Y, Tarr, l, dt, tend, str)
% 
% 
% 
% for t = 1:(tend/dt)-1
%               
%     VisTemperature(formfunction, X, Y, Tarr(:,:,t), l, str)    
% 
% end


function VisTemperature_transient(formfunction, X, Y, Tarr, l, str)


filename = 'transientconduction_surface.gif';
figure4 = figure;
figure(1)
set(gca, 'FontSize', 18)
set(gcf,'position', [34   164   560   420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
zlim([0 100]);


for i = 0:size(Tarr,3)-1
surf(X, Y, Tarr(:,:,i+1));

hold on
campos([58.6464  -61.7142  130.2245])
surf(X, -Y,Tarr(:,:,i+1));

zlabel('T')
xlabel('x')
ylabel('y')
zlim([0 100])
if nargin == 6
title(str)
end

    % Time gif
    drawnow
    frame = getframe(gcf);
    im=frame2im(frame);
    [imind,cm] = rgb2ind(im,64);
    if i == 0
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','writeMode','append','Delaytime',0)
    end
    hold off
end

filename2 = 'transientconduction_contour.gif';
figure(2)
set(gca, 'FontSize', 18)
set(gcf,'position', [613   163   560   420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
 set(gcf, 'color', 'white')
hold off

for i = 0:size(Tarr,3)-1
plot([0 linspace(0, l, 50) l], [0 formfunction(linspace(0, 1, 50)) 0])
hold on
plot([0 linspace(0, l, 50) l], [0 -formfunction(linspace(0, 1, 50)) 0])
contour(X,Y,Tarr(:,:,i+1))
contour(X,-Y,Tarr(:,:,i+1))
% contour(X,-Y,T)

ylim([-5, 5])
axis equal

    xlabel('x')
ylabel('y')
if nargin == 6
title(str)
end
colorbar
caxis([0 100])
drawnow
    frame = getframe(gcf);
    im=frame2im(frame);
    [imind,cm] = rgb2ind(im,64);
    if i == 0
        imwrite(imind,cm,filename2,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename2,'gif','writeMode','append','Delaytime',0)
    end
    hold off
end

filename3 = 'transientconduction_colour.gif';
figure(3)
set(gca, 'FontSize', 18)
set(gcf,'position', [613   163   560   420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
 set(gcf, 'color', 'white')
hold off


for i = 0:size(Tarr,3)-1

figure(4)
set(gca, 'FontSize', 18)
set(gcf,'position', [ 1190         162         560         420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
 set(gcf, 'color', 'white')
hold off
plot([0 linspace(0, l, 50) l], [0 formfunction(linspace(0, 1, 50)) 0])
hold on
plot([0 linspace(0, l, 50) l], [0 -formfunction(linspace(0, 1, 50)) 0])
pcolor(X,Y,Tarr(:,:,i+1))
shading interp
pcolor(X,-Y,Tarr(:,:,i+1))
shading interp
% pcolor(X,-Y,T)
% shading interp

ylim([-5, 5])
axis equal

xlabel('x')
ylabel('y')

if nargin == 6
title(str)
end

colorbar
caxis([0 100])
% 
drawnow
    frame = getframe(gcf);
    im=frame2im(frame);
    [imind,cm] = rgb2ind(im,64);
    if i == 0
        imwrite(imind,cm,filename3,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename3,'gif','writeMode','append','Delaytime',0)
    end
    hold off
end