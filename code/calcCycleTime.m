% function [cycle_time]=calcCycleTime(vep_Fr)

vep_Fr=compiledData.vep_Fr;
Fs=compiledData.Fs;
dur_in_sec=size(vep_Fr,4)./Fs;
epoch=1+(0:0.4:dur_in_sec).*Fs;
X_vep=(1:1:size(vep_Fr,4))/Fs;

for x=1:size(vep_Fr,1)
    for y=1:size(vep_Fr,2)
        vep_fr=squeeze(vep_Fr(x,y,:,:));
        for z=1:length(epoch)-1
            vep_epoch=vep_fr(:,epoch(z):epoch(z+1));
        end
    end
end