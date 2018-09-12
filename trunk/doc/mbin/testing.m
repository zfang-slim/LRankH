%% Analytic solutions
%
% test against analytic solutions, see <analytic_test.m |analytic_test.m|>

tmp = dlmread([resultsdir '/testing/error_mod.dat']);
h   = tmp(:,1);
e   = tmp(:,2);

% plot error
figure;
loglog(h,e,'k',h,.3*h.^2,'k--');xlim([2,20])
xlabel('h');ylabel('error');

% plot solutions
[D1,o,d,n] = odnread([resultsdir '/testing/D1.odn']);
[D2,o,d,n] = odnread([resultsdir '/testing/D2.odn']);
[z,x] = odn2grid(o,d,n);

D1 = reshape(D1,n);
D2 = reshape(D2,n);

iz = floor(length(z)/2) + 1;

figure;
plot(x,real(D1(iz,:)),'k',x,real(D2(iz,:)),'r',x,imag(D1(iz,:)),'k--',x,imag(D2(iz,:)),'r--');
legend('Re(numerical)','Re(analytic)','Im(numerical)','Im(analytic)');
xlabel('x [m]');ylabel('z [m]');
title(['solutions at z = ' num2str(z(iz)) 'm']);


[G1,o,d,n] = odnread([resultsdir '/testing/G1.odn']);
[G2,o,d,n] = odnread([resultsdir '/testing/G2.odn']);
[z,x] = odn2grid(o,d,n);
G1 = reshape(G1,n);
G2 = reshape(G2,n);

figure;
imagesc(1e-3*x,1e-3*z,real(G1),[-10 10]);axis equal tight;colormap(seiscol);
xlabel('x [km]');ylabel('z [km]');

figure;
imagesc(1e-3*x,1e-3*z,real(G2),[-10 10]);axis equal tight;colormap(seiscol);
xlabel('x [km]');ylabel('z [km]');

%% Jacobian test
% 
% Accuracy of the Jacobian, see <jacobian_test.m |jacobian_test.m|>

tmp = dlmread([resultsdir '/testing/error_jac.dat']);
h   = tmp(:,1);
e   = tmp(:,2);

% plot error
figure;
loglog(h,e,'k',h,1e3*h.^2,'k--');
xlabel('h');ylabel('error');

%% Adjoint test
%
% Adjoint test of the Jacobian, , see <adjoint_test.m |adjoint_test.m|>

adjoint_table = dlmread([resultsdir '/testing/error_adj.dat']);

adjoint_table(:,1:3) 

%% Scalability
par_table = zeros(5,3);
par_table(:,1) = [1 2 4 8 16]';
par_table(:,2) = [2459 1224 615 313 188]';
par_table(:,3) = par_table(1,2)./(par_table(:,2).*par_table(:,1));

par_table
