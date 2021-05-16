function make_copy_for_fieldtrip(D,csdFlag)



if csdFlag
    
    S = [];
    S.D = D;
    S.outfile = ['csdnanart_' D.fname]; %append 'c' on start of filename for CSD transformation
    
    tmpD = spm_eeg_copy(S);   
    
else
    S = [];
    S.D = D;
    S.outfile = ['nanart_' D.fname]; %append 'c' on start of filename for CSD transformation
    
    tmpD = spm_eeg_copy(S);

end

bs = tmpD.badsamples(:,:,:);
for ch = 1:size(bs,1)
    tmpD(ch,find(bs(ch,:)),1) = nan;
end
tmpD.save; clear tmpD





end