function e = JacobianTest(fh,x,dx,h)
% Test Jacobian using Taylor series
%
% e(h) = F(x+h*dx) - F(x) 0 J*dx = O(h^2)
%
% use:
%   e = JacobianTest(fh,x,dx,h)
%
% input:
%   fh - function handle [y,J] = fh(x)
%   x  - input vector
%   dx - perturbation
%   h  - vector with stepsizes
%
% output 
%   e  - vector with errors
%
%

[D0, J0] = fh(x);
dD = J0*dx;
D0 = gather(D0);
dD = gather(dD);
e = 0*h;
for k = 1:length(h)
    D1 = fh(x+h(k)*dx);
    D1 = gather(D1);
    e(k) = norm(D1 - D0 -h(k)*dD);
end