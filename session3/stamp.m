function [stecil] = stamp(i, j, X, Y)
%stecil calculate the linear equation for node (i, j)

%  input:
%    i         node number in x direction
%    j         node number in y direction
%    X         x position of the nodes
%    Y         y position of the nodes
%    b         right-hand side value for node (i,j)
%    alpha     alpha
%    Tinf      Tinf for Robin BC
%    boundary  defines the boundary conditions
%    verbose   verbositiy level
%
%  output:
%    stecil     linear equation for node (i,j)
%    b         new right-hand side value for node (i,j)


% Init

n = size(X, 1);
m = size(X, 2);
stecil = zeros(1, n*m);
index=@(ii, jj) ii + (jj-1)*n;

% Determine the node positon
if i==1
    nodePosition = 'north';
elseif i==size(X,1)
    nodePosition = 'south';
elseif j==1
    nodePosition = 'west';
elseif j==size(X,2)
    nodePosition = 'east';
else
    nodePosition = 'inner Node';   
end

% Calculate the equation for the correct node position
switch nodePosition
    
    case 'inner Node'
        
        
% Nomenclature:
%
%    NW(i-1,j-1)   Nw -  N(i-1,j) -  Ne     NE(i-1,j+1)
%
%                 |                 |
%
%       nW - - - - nw ------ n ------ ne - - - nE
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%   W(i, j-1) - - w - - P (i,j) - - e - -  E (i,j+1)
%                 |                 |
%       |         |        |        |       |
%                 |                 |
%      sW - - - - sw ------ s ------ se - - - sE
%
%                 |                 |
%
%   SW(i+1,j-1)   Sw  -  S(i+1,j)  - Se      SE(i+1,j+1)
%
% Indexing of stecil: 

%    D_4 - D_1 - D2
%     |     |     | 
%    D_3 - D_0 - D3
%     |     |     | 
%    D_2 -  D1 - D4



        %y_w=(Y(i,j)+Y(i, j-1))/2;
        %y_Sw=(Y(i+1, j) + Y(i+1, j-1))/2;
        %y_Se=(Y(i+1, j+1) + Y(i+1, j))/2;
        %y_e=(Y(i, j+1) + Y(i, j))/2;
        
        %y_s=(Y(i,j)+Y(i+1,j))/2;
        %y_sE=(Y(i+1, j+1) + Y(i, j+1))/2;
        %y_nE=(Y(i-1, j+1) + Y(i, j+1))/2;
        %y_n=(Y(i, j) + Y(i-1, j))/2;
        
        %y_Ne = (Y(i-1,j ) + Y(i-1, j+1))/2;
        %y_Nw = (Y(i-1, j-1) + Y(i-1, j))/2;
        
        %y_nW = (Y(i-1, j-1) + Y(i, j-1))/2;
        %y_sW = (Y(i, j-1) + Y(i+1, j-1))/2;
        
        y_NW = Y(i-1,j-1);   x_NW = X(i-1,j-1); % checked
        y_N  = Y(i-1,j  );   x_N  = X(i-1,j  ); % checked
        y_NE = Y(i-1,j+1);   x_NE = X(i-1,j+1); % checked
        y_W  = Y(i  ,j-1);   x_W  = X(i  ,j-1); % checked
        y_P  = Y(i  ,j  );   x_P  = X(i  ,j  ); % checked
        y_E  = Y(i  ,j+1);   x_E  = X(i  ,j+1); % checked
        y_SW = Y(i+1,j-1);   x_SW = X(i+1,j-1); % checked
        y_S  = Y(i+1,j  );   x_S  = X(i+1,j  ); % checked
        y_SE = Y(i+1,j+1);   x_SE = X(i+1,j+1); % checked
     
        y_Nw = (y_NW + y_N)/2;  x_Nw = (x_NW + x_N)/2;  % checked 
        y_Ne = (y_NE + y_N)/2;  x_Ne = (x_NE + x_N)/2;  % checked
        y_w  = (y_P  + y_W)/2;  x_w  = (x_P + x_W)/2;   % checked
        y_e  = (y_P  + y_E)/2;  x_e  = (x_P + x_E)/2;   % checked
        y_Sw = (y_S + y_SW)/2;  x_Sw = (x_S + x_SW)/2;  % checked
        y_Se = (y_SE + y_S)/2;  x_Se = (x_SE + x_S)/2;  % checked
        y_nW = (y_NW + y_W)/2;   x_nW = (x_NW + x_W)/2; % checked
        y_sW = (y_W  + y_SW)/2;  x_sW = (x_W  + x_SW)/2;% checked
        y_n  = (y_P  + y_N)/2;   x_n =  (x_P + x_N)/2;  % checked
        y_s  = (y_P  + y_S)/2;   x_s =  (x_P + x_S)/2;  % checked
        y_nE = (y_NE + y_E)/2;   x_nE = (x_NE + x_E)/2; % checked
        y_sE = (y_SE + y_E)/2;   x_sE = (x_SE + x_E)/2; % checked
        
        y_se = (y_Se + y_e)/2;  x_se = (x_Se + x_e)/2;  % checked
        y_sw = (y_Sw + y_w)/2;  x_sw = (x_Sw + x_w)/2;  % checked
        y_nw = (y_Nw + y_w)/2;  x_nw = (x_Nw + x_w)/2;  % checked
        y_ne = (y_Ne + y_e)/2;  x_ne = (x_Ne + x_e)/2;  % checked

        % Around s 
        
        dy_Sw_Se = y_Se - y_Sw;  dx_Sw_Se = x_Se - x_Sw;
        dy_Se_e = y_e - y_Se;    dx_Se_e = x_e - x_Se;
        dy_e_w = y_w - y_e;      dx_e_w = x_w - x_e;
        dy_w_Sw = y_Sw - y_w;    dx_w_Sw = x_Sw - x_w;
  
        % Around e

        dy_s_sE = y_sE - y_s;    dx_s_sE = x_sE - x_s;
        dy_sE_nE = y_nE - y_sE;  dx_sE_nE = x_nE - x_sE;
        dy_nE_n = y_n - y_nE;    dx_nE_n = x_n - x_nE;
        dy_n_s = y_s - y_n;      dx_n_s = x_s - x_n;

        % Around n
        
        dy_w_e   = y_e - y_w;    dx_w_e   = x_e - x_w;
        dy_e_Ne  = y_Ne - y_e;   dx_e_Ne  = x_Ne - x_e;
        dy_Ne_Nw = y_Nw - y_Ne;  dx_Ne_Nw = x_Nw - x_Ne;
        dy_Nw_w  = y_w - y_Nw;   dx_Nw_w  = x_w - x_Nw;
        
        
        % Around w
        dy_sW_s = y_s - y_sW;    dx_sW_s = x_s - x_sW;
        dy_s_n  = y_n - y_s;     dx_s_n  = x_n - x_s;
        dy_n_nW = y_nW - y_n;    dx_n_nW = x_nW - x_n;
        dy_nW_sW = y_sW - y_nW;  dx_nW_sW = x_sW - x_nW;
        

        % Around P
        
        dy_sw_se = y_se - y_sw; dx_sw_se = x_se - x_sw;
        dy_se_ne = y_ne - y_se; dx_se_ne = x_ne - x_se;
        dy_ne_nw = y_nw - y_ne; dx_ne_nw = x_nw - x_ne;
        dy_nw_sw = y_sw - y_nw; dx_nw_sw = x_sw - x_nw;

% Areas
    
        
        S_P                 = abs( x_ne*y_se - y_ne*x_se...
            + x_se*y_sw - y_se*x_sw...
            + x_sw*y_nw - y_sw*x_nw...
            + x_nw*y_ne - y_nw*x_ne) / 2;
        S_s                 = abs( x_e*y_Se - y_e*x_Se...
            + x_Se*y_Sw - y_Se*x_Sw...
            + x_Sw*y_w - y_Sw*x_w...
            + x_w*y_e - y_w*x_e) / 2;
        S_e                 = abs( x_nE*y_sE - y_nE*x_sE...
            + x_sE*y_s - y_sE*x_s...
            + x_s*y_n - y_s*x_n...
            + x_n*y_nE - y_n*x_nE) / 2;
        S_n                 = abs( x_Ne*y_e - y_Ne*x_e...
            + x_e*y_w - y_e*x_w...
            + x_w*y_Nw - y_w*x_Nw...
            + x_Nw*y_Ne - y_Nw*x_Ne ) / 2;
        S_w                =  abs( x_n*y_s - y_n*x_s...
            + x_s*y_sW - y_s*x_sW...
            + x_sW*y_nW - y_sW*x_nW...
            + x_nW*y_n - y_nW*x_n) / 2;          
        
        
        
        %$$$$$$$$$$$$$$$$$$$$$$ Stecil $$$$$$$$$$$$$$$$$$$

        build_inner

        %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        % East
        stecil(index(i, j+1)) = D3;  
        % West
        stecil(index(i, j-1)) = D_3;
        % South
        stecil(index(i+1, j)) = D1;
        % North
        stecil(index(i-1, j)) = D_1;
        % NW
        stecil(index(i-1, j-1)) = D_4;
        % NE
        stecil(index(i-1, j+1)) = D2;
        % SW
        stecil(index(i+1, j-1)) = D_2;
        % SE
        stecil(index(i+1, j+1)) = D4;
        % P
        stecil(index(i, j)) = D0;
        
    case 'south'
        stecil = zeros(1,size(X,2)*size(X,1));
        stecil(1,index(i,j)) = 1;
    case 'north'
        stecil = zeros(1,size(X,2)*size(X,1));
        stecil(1,index(i,j)) = 1;
    case 'east'
        stecil = zeros(1,size(X,2)*size(X,1));
        stecil(1,index(i,j)) = 1;
    case 'west'
        stecil = zeros(1,size(X,2)*size(X,1));
        stecil(1,index(i,j)) = 1;
end



