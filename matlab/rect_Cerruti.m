%--------------------------------------------------------------------------
function I = rect_Cerruti (X,Y,a,b)
% Finds the integrals  I = [I1; I2; I3; I4] needed to compute
% the displacements at the surface of an elastic half-space when
% subjected to a static, tangential, rectangular load of unit intensity
% q = 1 = Px/(2*a)/(2*b), where Px is the total horizontal load
% X Y are the coordinates relative to the center
% a b are the half-sizes of the rectangular load
%
% Written by Eduardo Kausel
%            MIT, 77 Mass. Ave., Room 1-271, Cambridge, MA 02139
%            Version 1, August 15, 2012
%
% Note: the actual displacements are
%       ux = q/pi/G *((1-nu)*I1 + nu*I2)
%       uy = q/pi/G * nu*I3
%       uz = q/pi/G *(1-nu)*I4
% where q = actual load intensity, G = shear modulus, nu = Poisson's ratio


x = X+a; % coordinates relative to the left lower corner
y = Y+b;
x2 = x/2;
y2 = y/2;
if abs(X)<=a & abs(Y)<=b
  % Normal point within rectangle
  if abs(X)~=a & abs(Y)~=b
    % arbitrary point
    I1 = integ1(x2,y2)+integ1(a-x2,y2)+integ1(x2,b-y2)+integ1(a-x2,b-y2);
    I2 = integ2(x2,y2)+integ2(a-x2,y2)+integ2(x2,b-y2)+integ2(a-x2,b-y2);
    I3 = integ3(x2,y2)-integ3(a-x2,y2)-integ3(x2,b-y2)+integ3(a-x2,b-y2);
    I4 = integ4(x2,y2)-integ4(a-x2,y2)+integ4(x2,b-y2)-integ4(a-x2,b-y2);
  elseif abs(X)~=a & abs(Y)==b
    % horizontal edge, not a corner
    if y==0
      I1 =  integ1(x2,b)+integ1(a-x2,b);
      I2 =  integ2(x2,b)+integ2(a-x2,b);
      I3 = -integ3(x2,b)+integ3(a-x2,b);
      I4 =  integ4(x2,b)-integ4(a-x2,b);
    else
      I1 = integ1(x2,b)+integ1(a-x2,b);
      I2 = integ2(x2,b)+integ2(a-x2,b);
      I3 = integ3(x2,b)-integ3(a-x2,b);
      I4 = integ4(x2,b)-integ4(a-x2,b);
    end
  elseif abs(X)==a & abs(Y)~=b
    % vertical edge, not a corner
    if x==0
      I1 =  integ1(a,y2)+integ1(a,b-y2);
      I2 =  integ2(a,y2)+integ2(a,b-y2);
      I3 = -integ3(a,y2)+integ3(a,b-y2);
      I4 = -integ4(a,y2)-integ4(a,b-y2);
    else
      I1 = integ1(a,y2)+integ1(a,b-y2);
      I2 = integ2(a,y2)+integ2(a,b-y2);
      I3 = integ3(a,y2)-integ3(a,b-y2);
      I4 = integ4(a,y2)+integ4(a,b-y2);
    end
  else
    % It is one of the corners
    I1 =  integ1(a,b);
    I2 =  integ2(a,b);
    if x==0 & y==0
      % left lower corner
      I3 =  integ3(a,b);
      I4 = -integ4(a,b);
    elseif x==0 & y>0
      % left upper corner
      I3 = -integ3(a,b);
      I4 = -integ4(a,b);
    elseif x>0 & y==0
      % right lower corner
      I3 = -integ3(a,b);
      I4 = integ4(a,b);
    elseif x>0 & y>0
      % right upper corner
      I3 = integ3(a,b);
      I4 = integ4(a,b);
    end
  end
  I = [I1;I2;I3;I4];
else
  % Point is beyond the rectangle
  if X<-a
    X = -X;
    x = X+a;
    x2 = x/2;
    sgnx = -1;
  else
    sgnx = 1; 
  end
  if Y<-b
    Y = -Y;
    y = Y+b;
    y2 = y/2;
    sgny = -1;
  else
    sgny = 1;
  end
  if abs(X)>a & abs(Y)>b
    % point beyond both the horizontal and vertical edges
    J1 = rect_Cerruti(x2,y2,x2,y2);
    J2 = rect_Cerruti(x2-a,y2-b,x2-a,y2-b);
    J3 = rect_Cerruti(x2,y2-b,x2,y2-b);
    J4 = rect_Cerruti(x2-a,y2,x2-a,y2);
    I = J1+J2-J3-J4;
  elseif abs(X)>a & abs(Y)<=b
    % point beyond left or right edge
    J1 = rect_Cerruti(x2,Y,x2,b);
    J2 = rect_Cerruti(x2-a,Y,x2-a,b);
    I = J1-J2;
  elseif abs(X)<=a & abs(Y)>b
    % point beyond upper or lower edge
    J1 = rect_Cerruti(X,y2,a,y2);
    J3 = rect_Cerruti(X,y2-b,a,y2-b);
    I = J1-J3;
  end
  if sgnx<0, I(3:4) = -I(3:4); end
  if sgny<0, I(3)   = -I(3);   end
end
endfunction

function I1 = integ1(a, b)
I1 = b*asinh(a/b)+a*asinh(b/a);
endfunction

function I2 = integ2(a, b)
I2 = b*asinh(a/b);
endfunction

function I3 = integ3(a, b)
I3 = a+b-sqrt(a^2+b^2);
endfunction

function I4 = integ4(a, b)
I4 = a*atan(b/a)+0.5*b*log(1+(a/b)^2);
I4 = -I4;
endfunction
