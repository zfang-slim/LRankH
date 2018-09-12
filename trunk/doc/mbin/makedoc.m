setpaths;
opts = [];
opts.outputDir = [basedir '/html'];
opts.format    = 'html';


%% modeling
publish('modeling',opts);

savefig(1,[basedir '/tex/figs_source/source']);
latextable([error1' error2' error3'],'name',[basedir '/tex/tables/source.tex'],'Horiz',{'linear','quadratic','cubic'},'Vert',{'0','1','2','3','4','5'},'format','%1.6e','Hline',[1,NaN]);

%% testing
publish('testing',opts);

savefig([1 2],[basedir '/tex/figs_source/modeling']);
savefig([3 4],[basedir '/tex/figs_source/vz']);
savefig(5,[basedir '/tex/figs_source/error_jacobian']);

latextable(adjoint_table(:,1:3),'name',[basedir '/tex/tables/adjoint.tex'],'Horiz',{'$\Re\langle A\mbf{x},\mbf{y} \rangle$','$\Re\langle \mbf{x},A^{\dagger}\mbf{y} \rangle$','error'},'format','%1.6e','Hline',[1,NaN]);

latextable(par_table(:,1:3),'name',[basedir '/tex/tables/scale.tex'],'Horiz',{'\# of procs','time [s]','efficiency'},'format','%1.2f','Hline',[1,NaN]);


%% examples
publish('examples',opts);

savefig([1:6],[basedir '/tex/figs_source/fwi']);

savefig([7:10],[basedir '/tex/figs_source/rtm']);

savefig([11:17],[basedir '/tex/figs_source/tomo']);


%% index
html_file = publish('index',opts);
open(html_file);