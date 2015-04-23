function [fi, fix, fiy, fiz, vol] = baseFun(T,iel,X)
%baseFun Basis functions of Hex1 elements
% Computes base functions and derivatives of HEX element in one or several
% points
% [fi, fix, fiy, fiz, vol] = baseHexP1(T,iel,X)
% T is the HEX mesh class
% iel is the element index number
% X is a set of coordinates at which to evaluate the basefunctions.

iv = T.Connectivity(iel,:);
xc = T.XC(iv); yc = T.YC(iv); zc = T.ZC(iv);

% Scalar evaluation of basis functions in point x y c
A=[1,xc(1),yc(1),zc(1);...
   1,xc(2),yc(2),zc(2);...
   1,xc(3),yc(3),zc(3);...
   1,xc(4),yc(4),zc(4)];

fim = (A\eye(4));


np = size(X,1);
o1 = ones(1,np);
z1 = zeros(1,np);

x = X(:,1); y = X(:,2); z = X(:,3);
x = x(:)'; y = y(:)'; z = z(:)';
fi=fim'*[o1;x;y;z];
fix=fim'*[z1;o1;z1;z1];
fiy=fim'*[z1;z1;o1;z1];
fiz=fim'*[z1;z1;z1;o1];

%% Compute volume
% cross product between a and b
a = [xc(2)-xc(1),yc(2)-yc(1),zc(2)-zc(1)];
b = [xc(3)-xc(1),yc(3)-yc(1),zc(3)-zc(1)];
v = [a(2).*b(3)-a(3).*b(2); 
     a(3).*b(1)-a(1).*b(3); 
     a(1).*b(2)-a(2).*b(1)].';
% scalar product between v and c;
c = [xc(4)-xc(1),yc(4)-yc(1),zc(4)-zc(1)].';

vol=abs(v*c)/6;

