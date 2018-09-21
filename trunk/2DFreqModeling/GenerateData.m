modeldir = '/math/home/fangzl/Model/overthrust3d';
datadir = '/scratch/Zhilong/Data/overthrust3d';
addpath('../');
addpath(modeldir);
addpath('/math/home/fangzl/Project/zfangtool');
Modelname = 'overthrust_vp.rsf';
curdir = pwd;
mfile = mfilename;
copyfile([mfile '.m'], datadir);
cd('../');
startup;



cd(modeldir);
[A o d n] = rsfread(Modelname);
A         = reshape(A,n);
cd(datadir);

model.o = [o(1) o(3)];
model.d = [d(1) d(3)];
model.n = [n(1) n(3)];
model.nb = [50 50];
model.freq = [2:17];
model.f0 = 10;
model.t0 = 0.01;
model.zsrc = 25;
model.xsrc = 25:25:(n(2)-2)*25;
model.zrec = 25;
model.xrec = 25:25:(n(2)-2)*25;
model.saveLU = 0;
model.saveWave = 0;
model.appWave  = 0;
model.uratio   = .1;

nsrc = length(model.xsrc);
nrec = length(model.xrec);
nfreq = length(model.freq);
od    = [model.xrec(1), model.xsrc(1), model.freq(1)];
dd    = [model.xrec(2)-model.xrec(1), model.xsrc(2)-model.xsrc(1), model.freq(2)-model.freq(1)];
nd    = [nrec,nsrc,nfreq];

for i = 172:floor(n(2)/2)
    fprintf('%03.0f\n',i)
    v1 = A(:,i,:);
    m1 = 1e6./v1(:).^2;
    Q  = eye(length(model.xsrc));
    D1 = F(m1,Q,model);
    D1     = gather(D1);
    D1     = reshape(D1, nrec, nsrc, nfreq);
    filename = sprintf('Data%03.0f.mat',i);
    WriteAllData(filename, D1, nd, dd, od);
end
cd(curdir)
