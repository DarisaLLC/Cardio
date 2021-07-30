function cardiomyocyte_3 (N, dL, L, W, t, Cs)
%**************************************************************************
%                 Mechanical model for a cardiomyocyte
% Arguments:
%   N = power index controlling the shape of the contact stresses p
%     > 0 -->  p = (x/a)^N           (e.g. if N=1 --> linear distribution) 
%     < 0 -->  p = (x/a)/sqrt(1-(x/a)^2) (Boussinesq-Cerruti distribution)
%           (observe that it produces a nearly linear displacement pattern)
%       Note: whatever the value of N, p is antisymmetric w.r.t the center
%   dL = cell elongation
%    L = cell length
%    W = cell width
%    t = cell thickness
%   Cs = shear wave velocity of gel
% The input parameters are assumed to be given in the CGS system (cm-g-s)
%
% Assumptions:
%   The cell is a thin rectangular elastic plate of negligible stiffness
%   The cell is bonded onto the surface of a gel, which is idealized as an
%     elastic half-space with uniform properties
%   As the cell contracts, it exerts tangential shearing stresses onto
%     the gel which are symmetric with respect to the cell's midplane.
%   The shearing stresses are assumed to be distributed according to some
%     rule and scaled to produce the observed elongation. The rule
%     controlling the shape of the shearing forces is defined by the
%     parameter N as follows:
%     N >= 0  ==> shearing stresses distributed according to (x/a)^N
%        < 0  ==> distribution according to 1/sqrt(1-(x/a)^2)
%
%   Dynamic effects are neglected. Instead, a Cerruti-type solution
%   consisting of a set of narrow rectangular loads applied onto the
%   surface of an elastic half-space is used.
%
%
%           Written by Eduardo Kausel
%           MIT, Room 1-271, Cambridge, MA, kausel@mit.edu
%           V-2, March 6, 2013
%
%**************************************************************************

% Default data in the absence of arguments
if nargin <6, Cs = 200.;  end   % shear wave velocity of gel [cm/s]
if nargin <5, t = 0.0002; end   % thickness [cm]
if nargin <4, W = 0.0030; end   % width [cm]
if nargin <3, L = 0.0100; end   % total length of cell [cm]
if nargin <2, dL = L;     end   % total elongation, cm
if nargin <1, N = 1;      end   % parameter controlling shear distribution
                                % (default is linear)

% Basic data
rho = 1.08; % mass density of gel [g/cm^3]
nu = 0.5;   % Poisson's ratio of gel
nx = 1000;  % number of rectangular elements in length (x) direction

% Inferred data
dx = L/nx;  % length of each element [cm]
a = dx/2;   % half-length of element [cm]
b = W/2;    % half-width of element  [cm]
G_gel = rho*Cs^2;   % shear modulus of gel [dynes/cm^2]
A_cell = W*L;       % area of cell
A = 4*a*b;          % area of rectangular element

% Assemble flexibility matrix F for the contact surface with the gel
nx1 = nx-1;
x = [0:dx:nx1*dx];  % x-coord center of nodes at left [cm]
xn = x+a;           % coord. of center of nodes
C = 1/(pi*G_gel*A); % factor 1/A --> unit force on cell --> intensity of q
F = zeros(nx,nx);
% Assemble first column
for i=1:nx
  I = C*rect_Cerruti (x(i), 0, a, b);
  F(i,1) = (1-nu)*I(1)+nu*I(2);
end
% clone that first column
for i=1:nx1
  for j=1:nx-i
    F(i+j,j+1) = F(i,1);
    F(j+1,i+j) = F(i+j,j+1);
  end
end
F(1,2:nx) = F(2:nx,1);

% Define load vector and solve the system
q = applied_shearing_stress (N, nx);  % see function below
u = F*q;           % displacements due to normalized tangential forces
scale = dL/2/u(1); % scaling factor needed to match observed elongation
u = scale*u;       % displacements now match elongation
p = scale*q;       % at this point, q are actual forces acting on elements,
                   % each of which has an area A=W*dx
Q = sum(p(1:nx/2));     % total actual net shearing force (dipole force)
M = -xn*p;              % moment of dipole
d = M/Q/L;              % moment arm of dipole, as fraction of cell length
S = Q/(W*L/2)*10^(-6);  % average contacts stress, in megadynes/cm^2
R = pi/4*G_gel*L*dL;    % force by simple dipole formula

% Print data and results
fprintf ('Length of cardiomyocyte               = %12.5f [cm]\n', L)
fprintf ('Width   "       "                     = %12.5f [cm]\n', W)
fprintf ('Thicknesss      "                     = %12.5f [cm]\n', t)
fprintf ('Observed elongation                   = %12.5f [cm]\n', dL)
fprintf ('Average contractive strain  e0=dL/L   = %12.5f\n', dL/L)
fprintf ('Total force by simple dipole formula  = %12.5f [dyne]\n', R)
fprintf ('Total reaction force exerted by gel   = %12.5f [dyne]\n', Q)
fprintf ('Average contact stress on gel         = %12.5f [M-dyne/cm^2]\n', S)
fprintf ('Total dipole moment                   = %12.5f [dyne-cm]\n', M)
fprintf ('Moment arm as fraction of cell length = %12.5f \n', d)

% Plot displacements along contact area
plot(xn*10^4,u*10^4)
grid on
title('Displacements along contact area [\mu]')
titx = 'Position of node [\mu]';
xlabel(titx)
pause

% Strains
e = diff(u)/dx;
plot(x(2:nx)*10^4,e)
title('Net axial strain in cell')
xlabel(titx)
grid on
pause

% Scaled contact stresses
plot(xn*10^4,p/A*10^(-6))
grid on
title('Shearing stresses acting on gel [M-dyne/cm^2]')
titx = 'Position of node [\mu]';
xlabel(titx)
pause

% Stresses in cell
s = cumsum(p)/(W*t);
plot(xn*10^4,s*10^(-6))
title('Net tensile stress within cell [M-dyne/cm^2]')
xlabel(titx)
grid on
pause
close
return

%--------------------------------------------------------------------------
function [q] = applied_shearing_stress(N, nx)
% Defines the shape of the normalized contact stresses being applied by the
% cell to the gel
if N==0
  % simple dipole
  q = zeros(nx,1);
  q(1)  =  1; % pair    ----> (left end)
  q(nx) = -1; %         <---- (right end)
  fprintf ('Shearing force is a simple dipole\n');
elseif N>0
  % polynomial distribution (linear, quadratic, whatever)
  n2 = nx/2;
  q = [0.5:n2-0.5]*2/nx;
  q = q.^N;
  q = [q(n2:-1:1),-q]';
  fprintf ('Shearing force changes as power N = %5.2f\n', N);
else
  % Distribution as in Boussinesq-Cerruti plate problem
  n2 = nx/2;
  q = [1:n2]-0.5;
  q = q*2/nx;
  q = q./sqrt(1-q.^2);
  q = [q(n2:-1:1),-q]';  
  fprintf ('Shearing force follows the Boussinesq shape\n');
end

%--------------------------------------------------------------------------
function [I] = rect_Cerruti (X,Y,a,b)
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
return

function [I1] = integ1(a, b)
I1 = b*asinh(a/b)+a*asinh(b/a);
return

function [I2] = integ2(a, b)
I2 = b*asinh(a/b);
return

function [I3] = integ3(a, b)
I3 = a+b-sqrt(a^2+b^2);
return

function [I4] = integ4(a, b)
I4 = a*atan(b/a)+0.5*b*log(1+(a/b)^2);
I4 = -I4;
return