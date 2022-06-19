function VisTemperature_transient(formfunction, X, Y, Tarr, l, h1, dt, str)

% make finer grid for interpolation
%x = 0:0.25:l;
%y = 0:0.25:h1;
% new grid for interpolation
%[Xq,Yq] = meshgrid(x,y);

filename = 'transientconduction_surface.gif';
figure(1)
set(gca, 'FontSize', 18)
set(gcf,'position', [34   164   560   420]);
set(gcf,'paperUnits','centimeters','paperPosition',[0 0 10 8],'paperSize',[10 8])
zlim([-100 100]);


for i = 0:size(Tarr,3)-1
surf(X, Y, Tarr(:,:,i+1));
%colormap('jet');

% interpolate values on finer grid
% Vq = interp2(X,Y,Tarr(:,:,i+1),Xq,Yq,'cubic');

% size(Vq)
% size(Xq)
% size(Yq)
% interpolate values

%surf(Xq,Yq,Vq);
surf(X, -Y, Tarr(:,:,i+1));
%colormap('jet');
hold on
campos([50  -50  100])

%surf(Xq, -Yq,Vq);

zlabel('p')
xlabel('x')
ylabel('y')
zlim([-100 100])
title([num2str(i*dt) ,' s'])
if nargin == 6

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
plot([0 linspace(0, l, 50) l], [0 formfunction(linspace(0, 1, 50)) 0],  'k-','LineWidth',1)
hold on
plot([0 linspace(0, l, 50) l], [0 -formfunction(linspace(0, 1, 50)) 0],  'k-','LineWidth',1)
contour(X,Y,Tarr(:,:,i+1))
%colormap('gray');
contour(X,-Y,Tarr(:,:,i+1))
%colormap('gray');
%contour(X,-Y,T)

ylim([-5, 5])
axis equal

xlabel('x')
ylabel('y')
title([num2str(i*dt) ,' s'])
if nargin == 6
title(str)
end
colorbar
caxis([-100 100])
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
pcolor(X,Y,Tarr(:,:,i+1))
colormap('gray')
hold on
shading interp
pcolor(X,-Y,Tarr(:,:,i+1))
colormap('gray')
%shading interp
plot([0 linspace(0, l, 50) l], [0 formfunction(linspace(0, 1, 50)) 0], 'k-','LineWidth',1)
% hold on
plot([0 linspace(0, l, 50) l], [0 -formfunction(linspace(0, 1, 50)) 0], 'k-','LineWidth',1)
% pcolor(X,Y,Tarr(:,:,i+1))
% %colormap('autumn');
% shading interp
% pcolor(X,-Y,Tarr(:,:,i+1))
%colormap('autumn');
shading interp
%pcolor(X,-Y,T)
% shading interp

ylim([-5, 5])
axis equal

xlabel('x')
ylabel('y')
title([num2str(i*dt) ,' s'])

if nargin == 6
title(str)
end

colorbar
caxis([-100 100])
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