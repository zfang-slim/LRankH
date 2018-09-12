%% Examples
%

%% FWI
% see <example_fwi.m example_fwi.m>
[v,o,d,n] = odnread([basedir '/data/bg_vp.rsf']);
[v0,o,d,n] = odnread([basedir '/data/bg_v0.rsf']);

[m1,o,d,n] = odnread([basedir '/results/fwi/m_ls_1.odn']);
[m2,o,d,n] = odnread([basedir '/results/fwi/m_ls_2.odn']);
[m3,o,d,n] = odnread([basedir '/results/fwi/m_ls_3.odn']);
[m4,o,d,n] = odnread([basedir '/results/fwi/m_ls_4.odn']);

[z,x] = odn2grid(o,d*1e-3,n);

figure;
imagesc(x,z,1e-3*reshape(v,n),[1.5 4.5]);
xlabel('x [km]');ylabel('z [km]');axis equal tight;

figure;
imagesc(x,z,1e-3*reshape(v0,n),[1.5 4.5]);
xlabel('x [km]');ylabel('z [km]');axis equal tight;

figure;
imagesc(x,z,reshape(real(1./sqrt(m1)),n),[1.5 4.5]);
xlabel('x [km]');ylabel('z [km]');axis equal tight;

figure;
imagesc(x,z,reshape(real(1./sqrt(m2)),n),[1.5 4.5]);
xlabel('x [km]');ylabel('z [km]');axis equal tight;

figure;
imagesc(x,z,reshape(real(1./sqrt(m3)),n),[1.5 4.5]);
xlabel('x [km]');ylabel('z [km]');axis equal tight;

figure;
imagesc(x,z,reshape(real(1./sqrt(m4)),n),[1.5 4.5]);
xlabel('x [km]');ylabel('z [km]');axis equal tight;

%% RTM
% see <example_rtm.m example_rtm.m>

[v,o,d,n]  = odnread([basedir '/data/marmvp.rsf']);
v = reshape(v,n);
v = [1500*ones(60,501); v(1:301,500 + [1:501])];
n = size(v);
S = opKron(opSmooth(n(2),200),opSmooth(n(1),200));

% smooth model
m  = 1e6./v(:).^2;
m0 = S*m;
dm = m - m0;
dm = reshape(dm,n);
v0 = reshape(1e3./sqrt(m0),n);
[dmt,o,d,n] = odnread([basedir '/results/rtm/dmt.rsf']);dmt=reshape(dmt,n);

[z,x] = odn2grid(o,d,n);

figure;
imagesc(x,z,v,[1500 4700]);colorbar;
xlabel('x [m]');ylabel('z [m]');
axis equal tight;

figure;
imagesc(x,z,v0,[1500 4700]);colorbar;
xlabel('x [m]');ylabel('z [m]');
axis equal tight;

figure;
imagesc(x,z,dm,[-1 1]*1e-1);colormap(gray);
xlabel('x [m]');ylabel('z [m]');
axis equal tight;

figure;
imagesc(x,z,dmt,[-1 1]*1e6);colormap(gray);
xlabel('x [m]');ylabel('z [m]');
axis equal tight;

%% tomo
%
[Dt,od,dd,nd] = odnread([basedir '/results/tomo/Dt.odn']);
Dt = reshape(Dt,nd);
[D0t,od,dd,nd] = odnread([basedir '/results/tomo/D0t.odn']);
D0t = reshape(D0t,nd);
[Ct,od,dd,nd] = odnread([basedir '/results/tomo/Ct.odn']);
Ct = reshape(Ct,nd);

T = odnread([basedir '/data/croswell_time.odn']);
T = reshape(T,nd(1:2));

[k1,o,d,n] = odnread([basedir '/results/tomo/k1.odn']);
k1 = reshape(k1,n);

[dm,o,d,n] = odnread([basedir '/data/crosswell_vp.odn']);
dm = reshape(dm,n);

[dmt,o,d,n] = odnread([basedir '/results/tomo/dmt.odn']);
dmt = reshape(dmt,n);

[zrec,zsrc,t] = odn2grid(od,dd,nd);
[z,x] = odn2grid(o,d,n);

figure;
imagesc(t,zrec,-squeeze(Dt(:,11,:)),[-1 1]*2e2);colormap(gray);
ylabel('z_r [m]');xlabel('t [s]');xlim([0 1]);
set(gca,'plotboxaspectratio',[1 2 1]);

figure;
imagesc(t,zrec,-squeeze(D0t(:,11,:)),[-1 1]*2e2);colormap(gray);
ylabel('z_r [m]');xlabel('t [s]');xlim([0 1]);
set(gca,'plotboxaspectratio',[1 2 1]);

figure;
imagesc(t-1,zrec,fftshift(squeeze(Ct(:,11,:)),2),[-1 1]*4e4);colormap(gray);
ylabel('z_r [m]');xlabel('\Delta t [s]');xlim([-.25 .25]);set(gca,'plotboxaspectratio',[1 2 1]);
tmp=fftshift(squeeze(Ct(:,11,:)),2);
[~,it] = max(tmp,[],2);
hold on;plot(t(it)-1,zrec,'r')

figure;
imagesc(zrec,zsrc,T);
xlabel('z_r [m]');ylabel('z_s [m]');colorbar

figure;
imagesc(x,z,k1);
xlabel('x [m]');ylabel('z [m]');

figure;
imagesc(x,z,dm,[1800 2500]);
xlabel('x [m]');ylabel('z [m]');colorbar;


figure;
imagesc(x,z,1e3./sqrt(.25+dmt),[1800 2500]);
xlabel('x [m]');ylabel('z [m]');

