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
