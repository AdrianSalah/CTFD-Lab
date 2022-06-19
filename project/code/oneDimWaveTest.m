
%% Plucked string example 
dx = .01 ; 
dt = .001 ; 

c = 5; 
L = 10 ; 
stopTime = 30 ; 

r = c*dt/dx ;
n = L/dx + 1 ;

% define initial conditions
current = 0.5-0.5*cos (2*pi/L*[0:dx:L]) ;
past = current ;

%plot ( [ 0:dx:L] , current )
%axis ( [ 0 L -2 2 ] )

for t=0:dt:stopTime
    % calculate future positions of string
    future(1) = 0;
    future(2:n-1) = r^2*( current (1:n-2)+current(3:n) ) + 2*(1-r^2)*current(2:n-1)-past(2:n-1);
    future(n) = 0;

    % update timestep
    past = current ;
    current = future ;

    % Plot the graph after every 10th frame
    if mod( t/dt , 10) == 0
        plot ( [0:dx:L] , current )
        axis ( [ 0 L -2 2 ] )
        pause ( .001 )
    end
end

%% Jump rope example

dx = .01 ; 
dt = .001 ; 
c = 1 ; 
L = 10 ; 
stopTime = 30 ; 

r = c*dt/dx ;
n = L/dx + 1 ;

current = zeros ( 1 , n ) ;
past = zeros ( 1 , n ) ;

for t=0: dt : stopTime
    
    future(1) = sin(2*t ) ;
    future(2:n-1) = r^2*( current ( 1:n-2)+current ( 3:n ) ) + 2*(1-r^2)*current(2:n-1) - past(2:n-1);
    future(n) = sin(3*t ) ;

    
    past = current ;
    current = future ;

     
     if mod( t /dt , 10) == 0
         plot ( [0:dx:L] , current )
         axis ( [0 L -8 8] )
         pause ( .001 )
     end
end
 
%% Reflected wave example


dx = .01 ; 
dt = .001 ; 
c = 1 ; 
L = 10 ; 
stopTime = 30 ; 

r = c *dt /dx ;
n = L/dx + 1 ;

 current = zeros ( 1 , n ) ;
 
 current (3/ dx : 4 / dx ) = -.5+.5* cos ( linspace (0 ,2*pi , 4 / dx-3/dx+1) ) ;
 past = current ;

 
 
 for t=0: dt : stopTime
     
     future ( 1 ) = 2* r^2* current ( 2 ) + 2*(1-r^2)* current ( 1 ) - past ( 1 ) ;
     future ( 2 : n-1) = r^2*( current ( 1 : n-2)+current ( 3 : n ) ) + 2*(1-r^2)* current ( 2 : n-1) - past ( 2 : n-1);
     future (n) = 2* r^2* current (n-1) + 2*(1-r^2)* current(n) - past (n ) ;

     
     past = current ;
     current = future ;

     
     if mod( t /dt , 10) == 0
         plot ( [ 0 : dx :L] , current )
         axis ( [ 0 L -1 1 ] )
         pause ( .001 )
     end
 end
