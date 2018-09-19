function x = GaussNewton(fh, x, opt)

    if isfield(opt,'write')
        write = opt.write;
    else
        write = 0;
    end

    if isfield(opt, 'NLitermax')
        NLitermax = opt.itermax;
    else
        NLitermax = 10;
    end

    if isfield(opt, 'Litermax')
        Litermax = opt.Litermax;
    else
        Litermax = 50;
    end


    for i = 1:NLitermax
        [f dD J] = fh(x);
        dx       = lsqrSOL(size(J,1), size(J,2), J, dD, 10^-10, [], [], [], Litermax, 1);
        x        = x - dx;
        fprintf('%03.0f,    %3.3e\n', i, f);
    end
