%% Teaser 1

% initialize
a = rand(10, 1);
b = rand(10, 1);

% long solution
% scalar product
s = 0;
for i = 1:length(a)
    s = s + a(i) * b(i);
end
% short solution

s2 = a'*b;

% Check your solution 
display(s - s2)


%% Teaser 2

% initialize
a = rand(10, 1);
b = rand(10, 1);

% long solution
c = zeros(size(a));
%elementwise multiplication
for i = 1:length(a)
    c(i) = a(i) * b(i);
end

% short solution
c2 = a.*b;

% Check your solution 
display(norm(c - c2))


%% Teaser 3

% initialize
a = rand(10, 1);
b = rand(10, 1);


% long solution
c = zeros(size(a));
% calculate elementwise product, b ist backwards array
for i = 1:length(a)
    c(i) = a(i) * b(end+1 -i);
    
end

% short solution
%c2 = a.*fliplr(b);

% Check your solution 
%display(norm(c - c2))

%% Teaser 4

% initialize
A = rand(10);

% long solution
vec = zeros(size(A, 1) * size(A, 2), 1);
% array holds columnwise matrix elements
for i = 1:size(A, 1)
    for j = 1:size(A, 2)
        vec(i + (j-1)*size(A, 2)) = A(i, j);
    end
end


%A = [1 2; 2 3]
%z = reshape(A, [length(A)^2, 1])


% short solution
%vec2 = reshape(A, [size(A,1)*size(A,2), 1]);
%or
vec2 = reshape(A, [length(A)^2, 1]);

% Check your solution 
display(norm(vec - vec2))


%% Teaser 5

B = zeros(10);

% initialize
b = rand(size(B, 1) * size(B, 2), 1);

% long solution

counter = 1;
for j = 1:size(B, 2) 
    for i = 1:size(B, 2)
        B(i,j) = b(counter);  
        counter = counter + 1;
    end
end

% short solution (2 lines, but no for-loop)
%B2 = ?;
%?;

% Check your solution 
%display(norm(B - B2))
