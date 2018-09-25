function x = GaussNewton(fh, x, opt)

    if isfield(opt,'write')
        write = opt.write;
    else
        write = 0;
    end

    if isfield(opt, 'NLitermax')
        NLitermax = opt.NLitermax;
    else
        NLitermax = 10;
    end

    if isfield(opt, 'Litermax')
        Litermax = opt.Litermax;
    else
        Litermax = 50;
    end

    if isfield(opt, 'xbound')
	xbound = opt.xbound;
    else
        xbound = [-10^10,10^10];
    end 


    for i = 1:NLitermax
        [f dD J] = fh(x);
        dx       = lsqrSOL(size(J,1), size(J,2), J, dD, 10^-10, [], [], [], Litermax, 1);
        x        = x - dx;
        idx      = find(x > xbound(2));
        x(idx)   = xbound(2);
        idx      = find(x < xbound(1));
        x(idx)   = xbound(1);
        fprintf('%03.0f,    %3.3e\n', i, f);
        save(['x_' num2str(i) '.mat'], 'x');
    end
