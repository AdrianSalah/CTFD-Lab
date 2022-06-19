function [X,Y,dummyNodes, IndexInletBoundary, IndexNorthBoundary] = setUpMesh(dimX, dimY, l, h, formfunction)

x = linspace(0,l,dimX);
y = linspace(h,0,dimY);
[X,Y] = meshgrid(x,y);

index = @(ii, jj) ii + (jj-1) * dimY;

dummyNodes = [];



% TODO: clean-up -> ugly implementation
IndexInletBoundary = 0;

for j=1:dimY
   if formfunction(0) > Y(j,1)
       IndexInletBoundary = IndexInletBoundary+1;
   end
end
IndexInletBoundary = dimY - IndexInletBoundary;

northBoundaryForm = double(formfunction(X(1,:)/l));
IndexNorthBoundary = zeros(1,dimX);

for i=1:dimX
   for j=1:dimY
       if Y(j,i) <= northBoundaryForm(i)
           IndexNorthBoundary(i) = index(j,i);
           break;
       end
   end
end


for col=1:dimX
    dummyNodes = [dummyNodes; double(Y(:,col) > formfunction(x(col)/l))];
end

end

