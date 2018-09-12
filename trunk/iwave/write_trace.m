function [] = write_trace(fname,Xs,Ys,Xr,Yr,Zs,Zr,T,data)
% Z = [oz dz nz];

segyid = fopen(fname,'w','b');


SuHeader.SegyFormatRevisionNumber=100;
SuHeader.DataSampleFormat=5;
SuHeader.Rev=GetSegyHeaderBasics;
SuHeader.ns=T(3);
SuHeader.nsOrig=T(3);
SuHeader.dt=T(2)*1e6;
SuHeader.dtOrig=T(2)*1e6;
SuHeader.FixedLengthTraceFlag=1;
SuHeader.NumberOfExtTextualHeaders=0;
it = 1;
for ks = 1:Ys(3)
    for ls = 1:Xs(3);
        for kr = 1:Yr(3)
            for lr = 1:Xr(3)
                SuTraceHeaders = InitSegyTraceHeader(T(3),1e6*T(1));
                
                SuTraceHeaders.SourceSurfaceElevation = -Zs;
                SuTraceHeaders.SourceX = Xs(1) + (ls-1)*Xs(2);
                SuTraceHeaders.SourceY = Ys(1) + (ks-1)*Ys(2);
                
                SuTraceHeaders.ReceiverGroupElevation = -Zr;
                SuTraceHeaders.GroupX = Xr(1) + (lr-1)*Xr(2);
                SuTraceHeaders.GroupY = Yr(1) + (kr-1)*Yr(2);
                if length(data)
                    PutSegyTrace(segyid,data(:,it),SuTraceHeaders,SuHeader);
                else
                    PutSegyTrace(segyid,zeros(T(3),1),SuTraceHeaders,SuHeader);
                end
                it = it + 1;
            end
        end
    end
end
fclose(segyid);
