setpaths;

% read model
[v,o,d,n]  = odnread([basedir '/data/marmvp.rsf']);

% model grid
model.o = o;
model.d = d;
model.n = n;

% absorbing boundary
model.nb = [50 50];

% source/receiver grid
model.zsrc = 25;
model.xsrc = linspace(0,11000,10);
model.zrec = 15;
model.xrec = 0:100:11000;

% frequencies
model.freq = 4.25:.25:20;

% wavelet
model.f0 = 10;
model.t0 = 0;

% source matrix
Q = speye(length(model.xsrc));

spmd,
    np = numlabs;
end
np = np{1};

fid = fopen('testing_parallel.log','r');
tmp = fgetl(fid);
fclose(fid);
k = strfind(tmp,'torque');
conf = tmp(k:end);

for k = 1:5
    tic;
    D = F(1e6./v.^2,Q,model);
    T(k)=toc;    
end


fid = fopen([resultsdir '/testing/time.dat'],'a');

fprintf(fid,'%s:%2d,%1.1f (%1.1f)\n',conf,np,mean(T),std(T));

fclose(fid);