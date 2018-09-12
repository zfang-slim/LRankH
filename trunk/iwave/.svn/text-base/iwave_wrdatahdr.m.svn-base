function [] = iwave_wrdatahdr(model,label)
% iwave_wrdatahdr(model,label)
%% source/receiver/time
if ~isfield(model,'ysrc')
    model.ysrc = 0;
end

if ~isfield(model,'yrec')
    model.yrec = 0;
end

ot = model.t(1);
dt = model.t(2) - model.t(1);
nt = length(model.t);

[orec,drec,nrec] = grid2odn(model.xrec,model.yrec,model.zrec);
[osrc,dsrc,nsrc] = grid2odn(model.xsrc,model.ysrc,model.zsrc);
%% headers

[Gxr,Gyr,Gxs,Gys] = ndgrid(model.xrec,model.yrec,model.xsrc,model.ysrc);

fid = fopen([label '.hdr'],'w');
fprintf(fid,'%f %f %f %f %f %f\n',[Gys(:)';Gxs(:)';Gyr(:)';Gxr(:)']);
fclose(fid);

cmd1 = sprintf('a2b < %s.hdr n1 = 6 > %s_hdr.bin',label,label);
cmd2 = sprintf('sunull nt=%d ntr=%d dt=%f > %s',nt,nrec(1)*nrec(2)*nsrc(1)*nsrc(2),dt,[label '_null.su']);
cmd3 = sprintf('sushw < %s_null.su infile=%s_hdr.bin key=sy,sx,gy,gx | sushw key=selev a=%f | sushw key=gelev a=%f > %s_hdr.su',label,label,-osrc(3),-orec(3),label);
system(['module load CWP; ' cmd1 ';' cmd2 ';' cmd3]);
system(['rm ' label '_null.su ' label '.hdr ' label '_hdr.bin']);

