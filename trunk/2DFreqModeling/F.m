function [D,J] = F(m,Q,model)
% Frequency domain FD modeling operator
%
% use:
%   [D,J] = F(m,Q,model)
% input:
%   m                 - vector with gridded squared slowness in [km^2/s^2]
%   Q                 - source matrix. size(Q,1) must match source grid
%                       definition, size(Q,2) determines the number of
%                       sources, if size(Q,3)>1, it represents a
%                       frequency-dependent source and has to be
%                       distributed over the last dimension.
%   model.{o,d,n}     - physical grid: z = ox(1) + [0:nx(1)-1]*dx(1), etc.
%   model.nb          - number of points to add for absorbing boundary
%   model.freq        - frequencies
%   model.f0          - peak frequency of Ricker wavelet, 0 for no wavelet.
%   model.t0          - phase shift [s] of wavelet.
%   model.{zsrc,xsrc} - vectors describing source array
%   model.{zrec,xrec} - vectors describing receiver array.
%
% output:
%   D  - Data cube (nrec x nsrc x nfreq) as (distributed) vector. nsrc  = size(Q,2);
%                                                                 nrec  = length(zrec)*length(xrec)
%                                                                 nfreq = length(freq)
%   J  - Jacobian as pSPOT operator
%%
% *Example*
% Below example defines a grid [0,1000]^2 with 10 m gridspacing. The
% velocity is 2 km/s. The sources and receivers coincide at x = [0:10:1000]
% and z = 10m.
%
%%
% model.o = [0 0 0];
% model.d = [10 10 1];
% model.n = [101 101 1];
% model.nb = [10 10 0];
% model.freq = [5 10 15 20];
% model.f0 = 0;
% model.zsrc = 10;
% model.xsrc = 0:10:1000;
% model.zrec = 10;
% model.xrec = 0:10:1000;
% m = .25*ones(prod(model.n),1);
% Q = speye(length(model.xsrc));
% D = F(m,Q,model);
%
%%

% Author: Tristan van Leeuwen
%         Seismic Laboratory for Imaging and Modeling
%         Department of Earch & Ocean Sciences
%         The University of British Columbia
%
% Date: February, 2012
%
% You may use this code only under the conditions and terms of the
% license contained in the file LICENSE provided with this source
% code. If you do not agree to these terms you may not use this
% software.

if nargin < 4
    dogather = 0;
end

if nargout < 2
    saveLU = 0;
    saveWave = 0;
else
    if isfield(model, 'saveWave')
        saveWave = model.saveWave;
    else
        saveWave = 0;
    end

    if saveWave == 0
        saveLU = 1;
        if isfield(model, 'saveLU')
            saveLU = model.saveLU;
        end
    else
        saveLU = 0;
    end
end

if saveWave == 1
    appWave = 0;
    if isfield(model, 'appWave')
        appWave = model.appWave;
    end
    uratio = 0.1;
    if isfield(model, 'uratio')
        uratio = model.uratio;
    end
end



% comp. grid
ot = model.o-model.nb.*model.d;
dt = model.d;
nt = model.n+2*model.nb;
[zt,xt] = odn2grid(ot,dt,nt);

% data size
nsrc   = size(Q,2);
nrec   = length(model.zrec)*length(model.xrec);
nfreq  = length(model.freq);

% define wavelet
w = exp(1i*2*pi*model.freq*model.t0);
if model.f0
    % Ricker wavelet with peak-frequency model.f0
    w = (model.freq).^2.*exp(-(model.freq/model.f0).^2).*w;
end

% mapping from source/receiver/physical grid to comp. grid
Pr = opKron(opLInterp1D(xt,model.xrec),opLInterp1D(zt,model.zrec));
Ps = opKron(opLInterp1D(xt,model.xsrc),opLInterp1D(zt,model.zsrc));
Px = opKron(opExtension(model.n(2),model.nb(2)),opExtension(model.n(1),model.nb(1)));

% model parameter: slowness [s/m] on computational grid.
nu = 1e-3*Px*sqrt(m);

% distribute frequencies according to standard distribution
freq = distributed(model.freq);
w    = distributed(w);

% check source matrix input
if (size(Q,3)==1)&&(isdistributed(Q))
    Q = gather(Q);
end

Prt = (sparsedouble(Pr))';

spmd
    if isfield(model,'numProcess')
	LASTN = maxNumCompThreads(model.numProcess);
    end 
    codistr  = codistributor1d(2,[],[nsrc*nrec,nfreq]);
    freqloc  = getLocalPart(freq);
    wloc     = getLocalPart(w);
    nfreqloc = length(freqloc);
    Dloc     = zeros(nrec*nsrc,nfreqloc);

    if saveLU > 0
        IHAll = {};
    end

    if saveWave > 0;
        UAll = [];
        VAll = [];
    end

    if size(Q,3)==1
        for k = 1:nfreqloc
            Hk  = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
            IHk = opHInv(Hk);
            Uk  = IHk*(wloc(k)*(Ps'*Q));
            Dloc(:,k) = vec(Pr*Uk);

            if saveLU > 0
                IHAll{k} = IHk;
            end

            if saveWave > 0
                UAll{k} = Uk;
                VAll{k} = IHk' * Prt;
                if appWave == 1
                    UAll{k} = opUapp(UAll{k}, nt, uratio);
                    VAll{k} = opUapp(VAll{k}, nt, uratio);
                end
            end
        end
    else
        Qloc = getLocalPart(Q);
        for k = 1:nfreqloc
            Hk  = Helm2D(2*pi*freqloc(k)*nu,ot,dt,nt,model.nb);
            IHk = opHInv(Hk);
            Uk  = IHk*(wloc(k)*(Ps'*Qloc(:,:,k)));
            Dloc(:,k) = vec(Pr*Uk);

            if saveLU > 0
                IHAll{k} = IHk;
            end

            if saveWave > 0
                UAll{k} = Uk;
                VAll{k} = IHk' * Prt;
                if appWave == 1
                    UAll{k} = opUapp(UAll{k}, nt, uratio);
                    VAll{k} = opUapp(VAll{k}, nt, uratio);
                end
            end
        end
    end
    D = codistributed.build(Dloc,codistr,'noCommunication');
end

% vectorize output, gather if needed
D = vec(D);

% construct pSPOT operator
model.saveWave = saveWave;
model.saveLU   = saveLU;
auxvar = [];
if saveWave == 1
    auxvar.UAll = UAll;
    auxvar.VAll = VAll;
end
if saveLU == 1
    auxvar.IHAll = IHAll;
end
J = oppDF(m,Q,model,auxvar);
