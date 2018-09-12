function [] = iwave_wrmodel(m,model,label)
%iwave_wrmodel(m,model,label)


odnwrite([label '_vp.rsf'],1./sqrt(m(:)),model.o,model.d,model.n);
odnwrite([label '_rho.rsf'],0*m + 1e3,model.o,model.d,model.n);