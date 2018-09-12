function j = runexp(label,mode,conf,poolsize,varargin)
% run batch job and save job information to file '<label>.log'
% 
% use:
%   j = runexp(label,mode,conf,poolsize,varargin)
% 
% input:
%   label - name of script (without .m)
%   mode  - 'i' interactive, 'b' batch
%   conf  - configuration
%   poolsize - poolsize
%   varargin - additional options passed to batch command
%
switch mode
    case 'i'
        matlabpool('open',conf,poolsize);
        eval(label);
        matlabpool close
    case 'b'
        j = batch(label,'configuration',conf,'matlabpool',poolsize,varargin{:});
        savejobdata(j, [label '.log'],conf);
    otherwise
        display('unknown mode');
end

function savejobdata(j,fname,conf)

fid = fopen(fname,'w');

jobno = j.ID;
jobid = j.pGetJobSchedulerData.pbsJobIds{1};
time  = j.SubmitTime;
file  = j.fileDependencies{1};

fprintf(fid,'%s, %d, %s, %s, %s \n',time,jobno,jobid,file,conf);
fprintf(1  ,'%s, %d, %s, %s, %s \n',time,jobno,jobid,file,conf);

fclose(fid);