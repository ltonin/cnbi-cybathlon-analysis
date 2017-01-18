function ConfMat = cnbiproc_commacc(RealGT, RealComm, PadTypeId, UseInd)

for gt=1:length(PadTypeId)
    for fc=1:length(PadTypeId)
        ConfMat(gt,fc) = 100* sum(RealGT==PadTypeId(gt) & RealComm==(PadTypeId(fc)+hex2dec('6000')) & UseInd)/sum(RealGT==PadTypeId(gt) & UseInd);
    end
end