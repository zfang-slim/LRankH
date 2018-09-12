function e = DotTest(A,N)
% performs dottest on matrix
%
% e = <A*x,y>/<x,A'*y>
% 
% for random vectors x and y
%
% use:
%   e = DotTest(A,N)
%
% input:
%   A - matrix
%   N - number of times to test
%
% output
%   e - vector with results

[m,n] = size(A);

e = zeros(N,1);

for k = 1:N
    x = randn(n,1);
    y = A*randn(n,1);
    
    e(k) = gather((A*x)'*y/(x'*(A'*y)));
end