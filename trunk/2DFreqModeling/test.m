function test_suite = test
% unit tests for 2DFDFDModeling
%
% Tristan van Leeuwen, 2011
% tleeuwen@eos.ubc.ca
%
initTestSuite;

function setup
matlabpool open local 3

function teardown
matlabpool close

function testAnalytic1
% test against Analytic solution for constant medium
% added: June 14, 2011

model.o = [0 0];
model.d = [10 10];
model.n = [101 151];
model.nb = [50 50];
model.freq = [10];
model.f0 = 0;
model.t0 = 0;
model.zsrc = 500;
model.xsrc = 500;
model.zrec = 0:10:1000;
model.xrec = 0:10:1500;

Q  = 1./prod(model.d);
v0 = 2000;
m  = 1e6*ones(prod(model.n),1)./v0.^2;

D1 = G([v0;0],Q,model);
D2 = F(m,Q,model);
D2 = gather(D2);

assertVectorsAlmostEqual(D1,D2,'relative',0.05);

function testAnalytic2
% test against Analytic solution for linear medium
% added: June 14, 2011

model.o = [0 0];
model.d = [10 10];
model.n = [101 151];
model.nb = [50 50];
model.freq = [10];
model.f0 = 0;
model.t0 = 0;
model.zsrc = 500;
model.xsrc = 500;
model.zrec = 0:10:1000;
model.xrec = 0:10:1500;

Q  = 1./prod(model.d);
v0 = 2000;
alpha = 0.5;
z  = odn2grid(model.o,model.d,model.n);
zz = repmat(z',1,model.n(2));
m  = 1e6./(v0 + alpha*zz(:)).^2;

D1 = G([v0;alpha],Q,model);
D2 = F(m,Q,model);
D2 = gather(D2);

assertVectorsAlmostEqual(D1,D2,'relative',0.1);

function testSource
% test frequency-dependent source
% added: Feb. 1, 2012

model.o = [0 0];
model.d = [10 10];
model.n = [51 51];
model.nb = [10 10];
model.freq = [5 10];
model.f0 = 0;
model.t0 = 0;
model.zsrc = 10;
model.xsrc = 0:50:500;
model.zrec = 10;
model.xrec = 0:50:500;
nsrc = length(model.xsrc);
nrec = length(model.xrec);
nfreq = length(model.freq);

Q = speye(nsrc);
m = 1e6*ones(prod(model.n),1)./2000.^2;

D  = F(m,Q,model);
D  = invvec(D,[nrec nsrc nfreq]);
DD = F(m,D,model);
DD = invvec(DD,[nrec nsrc nfreq]);

D  = gather(D);
DD = gather(DD);

for k = 1:nfreq
    tmpa = vec(D(:,:,k)*D(:,:,k));
    tmpb = vec(DD(:,:,k));
    assertVectorsAlmostEqual(tmpa,tmpb,'relative',0.05);
end

function testJacobian1
% test Jacobian
% added: June 14, 2011

model.o = [0 0];
model.d = [10 10];
model.n = [101 101];
model.nb = [20 20];
model.freq = [10 15];
model.f0 = 10;
model.t0 = 0.01;
model.zsrc = 15;
model.xsrc = 0:100:1000;
model.zrec = 10;
model.xrec = 0:5:1000;

v0 = 2000;
m = 1e6/v0.^2*ones(prod(model.n),1);
dm = .01*randn(model.n); dm([1:20 end-20:end],:) = 0;
dm(:,[1:20 end-20:end]) = 0; dm = dm(:);

Q = speye(length(model.xsrc));

h  = 10.^[0:-1:-5];
fh = @(x)F(x,Q,model);
e  = JacobianTest(fh,m,dm,h);

assertVectorsAlmostEqual(-diff(log10(e)),2*ones(size(diff(e))),'relative',0.05);

function testJacobian2
% test Jacobian with extended source
% added: Feb 1, 2012

model.o = [0 0];
model.d = [10 10];
model.n = [101 101];
model.nb = [20 20];
model.freq = [10 15];
model.f0 = 10;
model.t0 = 0.01;
model.zsrc = 15;
model.xsrc = 0:50:1000;
model.zrec = 10;
model.xrec = 0:50:1000;

v0 = 2000;
m = 1e6/v0.^2*ones(prod(model.n),1);
dm = .01*randn(model.n); dm([1:20 end-20:end],:) = 0;
dm(:,[1:20 end-20:end]) = 0; dm = dm(:);

Q = distributed(randn(length(model.xsrc),length(model.xsrc),length(model.freq)));

h  = 10.^[0:-1:-5];
fh = @(x)F(x,Q,model);
e  = JacobianTest(fh,m,dm,h);

assertVectorsAlmostEqual(-diff(log10(e)),2*ones(size(diff(e))),'relative',0.05);

function testAdjoint1
% Adjoint test of Jacobian
% added: June 14, 2011

model.o = [0 0];
model.d = [10 10];
model.n = [51 51];
model.nb = [10 10];
model.freq = [10 15];
model.f0 = 10;
model.t0 = 0.01;
model.zsrc = 15;
model.xsrc = 0:10:500;
model.zrec = 10;
model.xrec = 0:5:500;

v0 = 2000;
m  = 1e6/v0.^2*ones(prod(model.n),1);
Q  = speye(length(model.xsrc));
J0 = oppDF(m,Q,model);

e = DotTest(J0,5);

assertElementsAlmostEqual(real(e),ones(size(e)));

function testAdjoint2
% Adjoint test of Jacobian with extended source
% added: Feb 1, 2012

model.o = [0 0];
model.d = [10 10];
model.n = [51 51];
model.nb = [10 10];
model.freq = [10 15];
model.f0 = 10;
model.t0 = 0.01;
model.zsrc = 15;
model.xsrc = 0:50:500;
model.zrec = 10;
model.xrec = 0:50:500;

v0 = 2000;
m  = 1e6/v0.^2*ones(prod(model.n),1);
Q  = distributed(randn(length(model.xsrc),length(model.xsrc),length(model.freq)));
J0 = oppDF(m,Q,model);

e  = DotTest(J0,5);

assertElementsAlmostEqual(real(e),ones(size(e)));


function testDerivative
% test derivative of Matrix w.r.t. wavenumber
% added June 14, 2011

o = [0 0];
d = [10 5];
n = [101 51];
nb = [10 10];

f = 10;
v = 2e3*ones(n);
k = 2*pi*f./v;
dk = 1e-4*randn(n); 
dk(1:nb(1),:) = 0; dk(end-nb(1):end,:) = 0;
dk(:,1:nb(2)) = 0; dk(:,end-nb(2):end) = 0;

[H1,dH1] = Helm2D(k(:),   o,d,n,nb);
[H2,dH2] = Helm2D(k(:)+dk(:),o,d,n,nb);


assertVectorsAlmostEqual(nonzeros(H1 + dH1*spdiags(dk(:),0,prod(n),prod(n))),nonzeros(H2),'relative',1e-3);
