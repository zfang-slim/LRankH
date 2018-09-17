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
model.saveLU = 0;
model.saveWave = 0;
model.appWave  = 0;
model.uratio   = .1;

v0 = 2000 * ones(model.n);
v1 = v0;
v1(26:30,:) = 2500;
m0 = 1e6./v0(:).^2;
m1 = 1e6./v1(:).^2;
dm = m1 - m0;
Q  = eye(length(model.xsrc));
[D1 J] = F(m1,Q,model);
[D0 J] = F(m0,Q,model);
dD     = D1 - D0;
pm     = J' * dD;
dD1    = J * dm;

dD = gather(dD);
dD = reshape(dD,length(model.xrec), length(model.xsrc), length(model.freq));

figure;imagesc(reshape(pm,model.n))
figure;imagesc(real(squeeze(dD(:,:,1))))

norm(vec(dD1(:,:,1)))

% dm2 = lsqr(J, dD1, 1e-6,1000);
tic
dm2 = lsqrSOL(size(J,1), size(J,2), J, dD1, 10^-16, [], [], [], 1000, 1);
figure;imagesc(reshape(dm2,model.n))
toc
