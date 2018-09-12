% directoy stuff
setpaths;
cd([basedir '/results/fwi']);

% model
[v,o,d,n]  = odnread([basedir '/data/bg_vp.rsf']);
[v0,o,d,n] = odnread([basedir '/data/bg_v0.rsf']);
mref       = 1e6./v.^2;
m0         = 1e6./v0.^2;

% model params
model.o    = o; model.d = d; model.n = n;
model.nb   = [20 20];
model.freq = 2.5:.5:20;
model.zsrc = 20;
model.xsrc = 0:100:7000;
model.zrec = 10;
model.xrec = 0:25:7000;
model.f0   = 15;
model.t0   = 0;

nsrc  = length(model.xsrc);
nrec  = length(model.xrec);
nfreq = length(model.freq);

Q = speye(nsrc);

% make/read data
if 0
    D = F(mref,Q,model);
    [od,dd,nd] = grid2odn(model.xrec,model.xsrc,model.freq);
    odnwrite([basedir '/data/bg_data.odn'],gather(D),od,dd,nd);
else
    D = odnread([basedir '/data/bg_data.odn']);
    D = distributed(D);
end
D = invvec(D,[nrec*nsrc nfreq]);

% inversion
If = {[1:3],[7:9],[19:21],[31:33]};

odnwrite('m_ls_0.odn',m0,o,d,n);

for k = 1:length(If)
    modelk      = model;
    modelk.freq = model.freq(If{k});
    Dk          = vec(D(:,If{k}));
    
    fh = @(x)f_ls(x,Q,Dk,modelk);
    
    options.fid     = fopen(['iter_ls_' num2str(k) '.log'],'w');
    options.itermax = 10;
    
    m0 = odnread(['m_ls_' num2str(k-1) '.odn']);
    
    mn = mylbfgs(fh,m0,options);
    
    odnwrite(['m_ls_' num2str(k) '.odn'],mn,o,d,n);   
end


