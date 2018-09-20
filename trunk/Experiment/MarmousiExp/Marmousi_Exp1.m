modeldir = '/math/home/fangzl/Model/Marmousi-2';
datadir  = '/math/home/fangzl/Data/Marmousi';
resdir   = '/math/home/fangzl/Result/Marmousi/Exp1';
addpath('../../');
addpath(modeldir);
addpath('/math/home/fangzl/Project/zfangtool');
Modelname = 'marmosi_vp_ini_smooth.mat';
Dataname  = 'Data2-17Hz.mat';
curdir = pwd;

if ~exist(resdir, 'dir')
    mkdir(resdir);
end

mfile = mfilename;
copyfile([mfile '.m'], datadir);
cd('../../');
startup;

Ifreq = {[1:4],[4:7],[7:10],[10:13],[13:16]};



cd(modeldir);
[A n d o] = ReadAllData(Modelname);
A         = reshape(A,n);
cd(datadir);
[D nd dd od] = ReadAllData(Dataname);
cd(resdir);

D = squeeze(D);
[model.zrec, model.xrec, model.zsrc, model.xsrc, model.freq] = odn2grid(od, dd, nd);

model.o = o;
model.d = d;
model.n = n;
model.nb = [50 50];
model.f0 = 10;
model.t0 = 0.01;
model.saveLU = 0;
model.saveWave = 1;
model.appWave  = 1;
model.uratio   = .01;

nsrc = length(model.xsrc);
nrec = length(model.xrec);
nfreq = length(model.freq);
od    = [model.xrec(1), model.xsrc(1), model.freq(1)];
dd    = [model.xrec(2)-model.xrec(1), model.xsrc(2)-model.xsrc(1), model.freq(2)-model.freq(1)];
nd    = [nrec,nsrc,nfreq];

v0 = A(:,:);
m0 = 1e6./v0(:).^2;
Q  = eye(length(model.xsrc));

opt.Write = 1;
opt.NLitermax = 2;
opt.Litermax  = 10;

for k = 1:length(Ifreq)
    modelk      = model;
    modelk.freq = model.freq(Ifreq{k});
    Dk          = D(:,:,Ifreq{k});
    Dk          = distributed(Dk);
    Dk          = vec(Dk);
    fh          = @(x) misfit_GN(x, Dk, Q, modelk);
    m0          = GaussNewton(fh, m0, opt);
    vk          = reshape(1./sqrt(m0),model.n);
    filename    = ['v_' num2str(k) '.mat'];
    WriteAllData(filename, vk, model.n, model.d, model.o);
    Expsub      = ['xm' num2str(k)];
    mkdir(Expsub);
    movefile('*.mat', Expsub);
end

filename    = ['vfinal.mat'];
WriteAllData(filename, vk, model.n, model.d, model.o);

cd(curdir)
