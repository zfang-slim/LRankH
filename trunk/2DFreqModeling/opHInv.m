classdef opHInv < opSpot
% SPOT operator: Inverse of the Helmholtz matrix with LU decomposition
%
% use:
%   IH = opHInv(A)
%
% input:
%   A - Helmholtz matrix
%
% output:
%   IH - Inverse of the Helmholtz matrix with LU factroziations
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
        LL, UU, Pp, Qp, R;
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Constructor
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function op = opHInv(A)

           m = size(A,1);
           n = size(A,2);

           op = op@opSpot('opHInv', m, n);
           op.cflag     = 1;
           op.linear    = 1;
           op.children  = [];
           op.sweepflag = true;
           op.m         = m;
           op.n         = n;
           op.LL        = [];
           op.UU        = [];
           op.Pp        = [];
           op.Qp        = [];
           op.R         = [];

           [op.LL, op.UU, op.Pp, op.Qp, op.R] = lu(A);

       end %constructor



    end


    methods ( Access = protected )
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Multiply
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       function out = multiply(op,in,mode)
           if mode==1
               out = op.Qp*(op.UU\(op.LL\(op.Pp*(op.R\(in)))));
           else
               out = op.R'\(op.Pp'*(op.LL'\(op.UU'\(op.Qp'*in))));
           end

       end %multiply

    end %protected methods

end %classdef
