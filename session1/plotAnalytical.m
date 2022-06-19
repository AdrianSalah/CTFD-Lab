function  plotAnalytical(resultArray, xRange, t_end, delta_t, caseName)
% function plots analytical solution stored in resultArray and exports
% video in a gif-file
%   Parameters:
%
%   resultArray:    results stored row-wise for every timestep
%   xRange:         spatial points for which we solve the function
%   t_end:          final time 
%   delta_t:        timestepsize
%   caseName:       string describing the current case 

gcf = figure;
axis tight manual % this ensures that getframe() returns a consistent size
%filename = 'Convection.gif';
filename = strcat(caseName, '.gif');

for t = 0:(t_end/delta_t)
              
    plot(xRange,resultArray(t+1,:))
    ylim([-1 1])
    xlabel('x in m')
    ylabel('\theta ')    
    title([num2str(t*delta_t) ,' s'])
    
    drawnow 
      % Capture the plot as an image 
      frame = getframe(gcf); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,8); 
      % Write to the GIF File 
      if t == 0 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','Delaytime',0); 
      end 
end

end

