modeldir = '/math/home/fangzl/Model/Marmousi-2';
datadir = '/scratch/Zhilong/Data/Marmousi';
addpath('../');
addpath(modeldir);
addpath('/math/home/fangzl/Project/zfangtool');
Modelname = 'marmosi_vp.mat';
curdir = pwd;
mfile = mfilename;
copyfile([mfile '.m'], datadir);
cd('../');
startup;



cd(modeldir);
[A n d o] = ReadAllData(Modelname);
A         = reshape(A,n);
cd(datadir);

model.o = o;
model.d = d;
model.n = n;
model.nb = [50 50];
model.freq = [2:17];
model.f0 = 10;
model.t0 = 0.01;
model.zsrc = 20;
model.xsrc = 50:50:10950;
model.zrec = 20;
model.xrec = 50:50:10950;
model.saveLU = 0;
model.saveWave = 0;
model.appWave  = 0;
model.uratio   = .1;

nsrc = length(model.xsrc);
nrec = length(model.xrec);
nfreq = length(model.freq);
od    = [model.zrec(1), model.xrec(1), model.zsrc(1), model.xsrc(1), model.freq(1)];
dd    = [1, model.xrec(2)-model.xrec(1), 1, model.xsrc(2)-model.xsrc(1), model.freq(2)-model.freq(1)];
nd    = [1,nrec,1,nsrc,nfreq];


v1 = A(:,:);
m1 = 1e6./v1(:).^2;
Q  = eye(length(model.xsrc));
D1 = F(m1,Q,model);
D1     = gather(D1);
D1     = reshape(D1, [1, nrec, 1, nsrc, nfreq]);
filename = sprintf('Data2-17Hz.mat');
WriteAllData(filename, D1, nd, dd, od);

cd(curdir)
