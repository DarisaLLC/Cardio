function x = cardiomyocyte (N, dL, L, W, t, Cs)
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
dx = L/nx  % length of each element [cm]
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


return


