function test_suite = test
% test for other misc. functions
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
initTestSuite;

function testodn
% basic test of odn2grid <-> grid2odn pair
% added June 14, 2011

o = [0 0 0];
d = [10 20 30];
n = [10 5 3];

[z,x,y] = odn2grid(o,d,n);
[o1,d1,n1] = grid2odn(z,x,y);

assertEqual(o,o1);
assertEqual(d,d1);
assertEqual(n,n1);

z = 0:10:1000;
x = 20:5:500;
y = 500:pi:800;

[o,d,n] = grid2odn(z,x,y);
[z1,x1,y1] = odn2grid(o,d,n);

assertEqual(x,x1);
assertElementsAlmostEqual(y,y1);
assertEqual(z,z1);

function testodnio
% basic test of odnread <-> odnwrite pair
% added June 14, 2011

z = 0:10:1000;
x = 20:5:500;
y = 500:10:800;

[o,d,n] = grid2odn(z,x,y);

a = randn(n) + 1i*randn(n);

odnwrite('test.odn',a,o,d,n);

[a1,o1,d1,n1] = odnread('test.odn');
[z1,x1,y1] = odn2grid(o1,d1,n1);

assertElementsAlmostEqual(x,x1);
assertElementsAlmostEqual(y,y1);
assertElementsAlmostEqual(z,z1);
assertElementsAlmostEqual(a(:),a1(:),'relative',1e-6);

delete('test.odn*')

function testVec1
% basic test of distvec and distunvec
% added June 21, 2011

n = [20 5 11];

X = distributed.randn(n);
x = vec(X);
Y = invvec(x,n);

assertEqual(X,Y);

x = distributed.randn(prod(n),1);
X = invvec(x,n);
y = vec(X);

assertEqual(x,y);


function testVec2
% basic test of distvec and distunvec
% added June 21, 2011

n = [20 5 2];

X = distributed.randn(n);
x = vec(X);
Y = invvec(x,n);

assertEqual(X,Y);

x = distributed.randn(prod(n),1);
X = invvec(x,n);
y = vec(X);

assertEqual(x,y);

function testNorms2
% gradient tests of norm functions
% added June 14, 2011

x  = randn(100,1);
dx = randn(100,1);
dx = dx/norm(dx);

s = abs(randn(100,1));
J{1} = @(x)twonorms(x,s);
m = randi(10,1);
J{2} = @(x)hubers(x,m);
df = randi(10,1);
J{3} = @(x)students(x,df);
J{4} = @(x)hybrid(x,m);

h = 10.^[0:-1:-5];
for k = 1:length(J)
    e = GradTest(J{k},x,dx,h);
    
    assertVectorsAlmostEqual(-diff(log10(e)),2*ones(1,length(h)-1),'relative',1e-2);
end