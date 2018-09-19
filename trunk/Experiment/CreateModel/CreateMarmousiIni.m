modeldir  = '/math/home/fangzl/Model/Marmousi-2';
modelfile = 'marmosi_vp.mat';
modelfileout = 'marmosi_vp_ini.mat';
curdir    = pwd;
cd(modeldir);

[A n d o] = ReadAllData(modelfile);
smooth_level = 150;
S1 = opSmooth(n(1), smooth_level);
S2 = opSmooth(n(2), smooth_level);
Aout = S1 * A;
Aout = (S2 * Aout')';
avec = mean(Aout, 2);
Aout = repmat(avec, 1, n(2));


figure;imagesc(A)
figure;imagesc(Aout);
WriteAllData(modelfileout, Aout, n, d, o);
cd(curdir)
