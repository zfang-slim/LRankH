function e = GradTest(fh,x,dx,h)
% test gradient using Taylor series
%
% e*(h) = f(x + h*dx)- f(x) - h*g'*dx
%
% use
%   e = GradTest(fh,x,dx,h)
% 
% input:
%   fh - function handle [f,g] = fh(x)
%   x  - input vector
%   dx - perturbation
%   h  - vector with stepsizes
%
% ouput:
%   e  - vector with errors

[f,g] = fh(x);
e = 0*h;
for k = 1:length(h)
    e(k) = abs(fh(x+h(k)*dx) - f - h(k)*g'*dx);
end