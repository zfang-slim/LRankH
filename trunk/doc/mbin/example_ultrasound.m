% model
[x,o,d,n] = odnread('../data/shepp-logan.odn');

% parameters
ng = n(1);
nr = 20;
nt = 20;
f  = 20;
v0 = 2000;

% operator
K = opFunction(nr*nt,ng*ng,@(x,mode)pw_scat(x,mode,ng,nr,nt,v0,f),1,1);

I = speye(128);
R = kron(I(1:ng,:),I(1:ng,:));
W = opWavelet(128,128)*opMatrix(R');

% data
d = K*x;

% noise
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
d = d + .5*randn(size(d));

% ls reconstruction
xls = lsmr(K,d,[],[],[],[],10);

% l1 reconstruction
opts = spgSetParms('iterations',100,'bpTol',1e-1);
y    = spgl1(K*W',d,0,norm(K*xls-d),[],opts);
xl1  = W'*y;

%
mkdir('../results/ultrasound');
odnwrite('../results/ultrasound/xls.odn',xls,[0 0],1e-3*[1 1],ng*[1 1]);
odnwrite('../results/ultrasound/xl1.odn',xl1,[0 0],1e-3*[1 1],ng*[1 1]);


% plot
figure;
subplot(1,3,1);
imagesc(reshape(x,ng,ng),[0 1]);colormap(gray)
subplot(1,3,2);
imagesc(reshape(xls,ng,ng),[0 1]);colormap(gray)
subplot(1,3,3);
imagesc(reshape(xl1,ng,ng),[0 1]);colormap(gray)