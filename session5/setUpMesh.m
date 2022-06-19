function [X, Y] = setUpMesh(T, l, formfunction)
%T contains information about dimX and dimY
Y = zeros(size(T));


    x = linspace(0,l,size(T,2));
    xnorm = x/l;
    for i = 1:size(T,2)
        Y(:,i) = linspace(formfunction(xnorm(i)),0,size(T,1));

    end
    y = zeros(size(T,1),1);
    
    X= meshgrid(x,y);
    
   


    






end