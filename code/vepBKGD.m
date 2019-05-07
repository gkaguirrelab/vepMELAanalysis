function [ttf_bkgd]=vepBKGD(processedVEPdata,Fs,dur_in_sec)

        L=dur_in_sec*Fs;
        bkgd_freq=[2 5 11.25 22.5 45];
        
        background=cat(1,processedVEPdata(1).vep_bkgd,processedVEPdata(2).vep_bkgd,processedVEPdata(3).vep_bkgd);
        background2=cat(1,processedVEPdata(1).vep_bkgd,processedVEPdata(2).vep_bkgd,processedVEPdata(3).vep_bkgd);
        counter=0;
        for yyy=1:size(background,1)-1
            if isnan(background(yyy))==1
                background2=cat(1,background2(1:yyy-1-counter,:),background(yyy+1:end,:));
                counter=counter+1;
            end
        end
        
        for i=1:1000
            rand_trial=sort(randi(size(background2,1),1,21));
            background2=background2(rand_trial,:);
            backgroundM=nanmedian(background2,1);

            Y=fft(backgroundM);
            P2=abs(Y/L);
            P1=P2(:,1:L/2+1);
            P1(:,2:end-1)=2*P1(:,2:end-1);
            ttf_BKGD=P1;
            clear P1 P2 Y
        
            [fooof_bkgd]=runFOOOF_bkgd(ttf_BKGD);
            
            xdata=fooof_bkgd.freqs;
            ydata=10.^(fooof_bkgd.power_spectrum)-10.^(fooof_bkgd.bg_fit);

            for a=1:length(bkgd_freq)
                peak_freq=bkgd_freq(a);
                temp=abs(xdata-peak_freq);
                temp2=find(temp==min(temp));
                if length(temp2)>1
                    temp3=max(ydata(temp2));
                    peak_freq_loc(a)=find(ydata==temp3);
                else
                    peak_freq_loc(a)=temp2;
                end
               fooof_Fr(i,a)=ydata(peak_freq_loc(a));
            end
        end
        
        ttf_bkgd.bkgd_freq=bkgd_freq;
        ttf_bkgd.bkgd_fooof_fr=nanmedian(fooof_Fr,1);
end