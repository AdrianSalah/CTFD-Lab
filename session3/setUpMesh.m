function [X,Y] = setUpMesh(T, l, formfunction)

size_x = size(T,2);
size_y = size(T,1);

x = linspace(0,l,size_x);
X = meshgrid(x,linspace(0,1,size_y));

Y = zeros(size_y, size_x);
for i=1:1:size_x
    Y(:,i) = linspace(formfunction(x(i)/l),0,size_y);
end

end

