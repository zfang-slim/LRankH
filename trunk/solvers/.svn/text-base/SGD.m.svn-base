function [w] = SGD(fh,w,options)
% Simple incremental gradient method, basically a steepest descent method
% where a different (randomly chosen) term of the function is used at each iteration.
%
% use:
%   x = SGD(fh,x0,options)
% 
% input:
%   fh                  - function handle [f,g] = fh(x,i) where i is a scalar
%                         in [1,nInstances] indicating which term to evaluate
%   x0                  - initial guess
%   options.m           - total batchsize REQUIRED
%   options.itermax     - maximum number of full evaluations (100)
%   options.theta       - stepsize (1e-6)
%   options.fid         - fid for log file (1)
%   options.seed        - random seed.
%   options.iter0       - start of iteration count, usefull for restarting
%                          previous run (0).
%   options.write       - write iterates to disk (false)
%
% output:
%   x - final iterate

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

itermax = 100;
theta   = 1e-6;
fid     = 1;
seed    = 1;
iter0   = 0;
write   = 0;

if isfield(options,'m')
    m = options.m;
else
    fprintf(fid,'max batchsize options.m is a required parameter!\n');
    return;
end
if isfield(options,'itermax')
    itermax    = options.itermax;
end
if isfield(options,'theta')
    theta      = options.theta;
end
if isfield(options,'fid')
    fid        = options.fid;
end
if isfield(options,'write')
    write      = options.write;
end
if isfield(options,'iter0')
    iter0      = options.iter0;
end
if isfield(options,'seed')
    seed       = options.seed;
end
maxSGDiter = m*itermax;
stepSize   = theta;

s = RandStream.create('mt19937ar','seed',seed);
iter = iter0;
while iter < maxSGDiter

    i = randi(s,m,1);
    [f,g] = fh(w,i); 
 
    w = w - stepSize*g;
    
    iter = iter + 1;
    if mod(iter,m) == 0
        fprintf(fid,'%5d, %5d, %5d, %10.5e, %10.5e, %10.5e\n',iter,iter,1,stepSize,f,norm(g));
        if write
            dlmwrite(['x_' num2str(iter) '.dat'],w);
        end
    end
    
end
