classdef opUapp < opSpot
% SPOT operator: Approximate wavefield
%
% use:
%   opU = opUapp(U, ratio)
%
% input:
%   U - Collection of all wavefields, each colume represents a wavefield
%
% output:
%   opU - Low rank representation of all the wavfield
%

% Author: Zhilong Fang
%         IMAGING AND COMPUTING GROUP
%         Department of Mathematics
%         Earth Resources Laboratory
%         Massachusetts Institute of Technology
%
% Date: Septermber, 2018
%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        UU, SS, VV
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opUapp(U, size_u, ratio)
           
           if nargin < 3
               ratio = 0.1;
           end
           
           if issparse(U)
               U = full(U);
           end

           m = size(U,1);
           n = size(U,2);

           op = op@opSpot('opUapp', m, n);
           op.cflag     = 1;
           op.linear    = 1;
           op.children  = [];
           op.sweepflag = true;
           op.m         = m;
           op.n         = n;
           op.SS        = [];
           op.UU        = [];
           op.VV        = [];
           n_rank       = floor(min(size_u) * ratio);
           
           

           for i = 1:size(U,2)
               u            = reshape(U(:,i), size_u);
               [Ut St Vt]   = svd(u); 
               op.UU(:,:,i) = Ut(:,1:n_rank);
%                St           = diag(St);
               op.SS(:,:,i)   = St(1:n_rank,1:n_rank);
               op.VV(:,:,i)   = Vt(:,1:n_rank)';
           end

       end %constructor
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Elementwise product
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function out = times(op,in)
           n_in = size(in);
           if n_in(1) ~= op.m || n_in(2) ~= op.n
               error('The matrix dimension does not match.');
           else
               for i = 1 : op.n
                   u = vec(op.UU(:,:,i) * op.SS(:,:,i) * op.VV(:,:,i));
                   in(:,i) = u .* in(:,i);
               end            
           end
           out = in;    

       end %elementwise multiply



    end


    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply is not implemented 
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       function out = multiply(op,in,mode)
           if mode==1
               out = zeros(op.m, size(in,2));
               for i = 1:op.n
                   u   = vec(op.UU(:,:,i) * op.SS(:,:,i) * op.VV(:,:,i));
                   out = out + u * in(i,:);
               end
           else
               out = zeros(op.n, size(in,2));
               for i = 1:op.n
                   u   = vec(op.UU(:,:,i) * op.SS(:,:,i) * op.VV(:,:,i));
                   out(i,:) = u'* in;
               end
           end

       end %multiply
       
       
      

    end %protected methods

end %classdef
