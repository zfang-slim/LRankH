modeldir      = '/math/home/fangzl/Model/Marmousi-2';
modelfile     = 'marmosi_vp.mat';
modelfileout1 = 'marmosi_vp_ini_smooth.mat';
modelfileout2 = 'marmosi_vp_ini_1D.mat';
curdir        = pwd;
cd(modeldir);

[A n d o]    = ReadAllData(modelfile);
smooth_level = 150;
S1 = opSmooth(n(1), smooth_level);
S2 = opSmooth(n(2), smooth_level);
Aout  = S1 * A;
Aout  = (S2 * Aout')';
Aout1 = Aout;
avec  = mean(Aout, 2);
Aout2 = repmat(avec, 1, n(2));


figure;imagesc(A)
figure;imagesc(Aout1);
figure;imagesc(Aout2);
WriteAllData(modelfileout1, Aout1, n, d, o);
WriteAllData(modelfileout2, Aout2, n, d, o);
cd(curdir)
