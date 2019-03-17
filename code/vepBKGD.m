function [ttf_bkgd_Fr,ttf_bkgdCI_Fr,vep_bkgd]=vepBKGD(processedVEPdata,Fs,A,norm_vep)

        background=cat(1,processedVEPdata(1).vep_bkgd,processedVEPdata(2).vep_bkgd,processedVEPdata(3).vep_bkgd);
        background2=cat(1,processedVEPdata(1).vep_bkgd,processedVEPdata(2).vep_bkgd,processedVEPdata(3).vep_bkgd);
        counter=0;
        for yyy=1:size(background,1)-1
            if isnan(background(yyy))==1
                background2=cat(1,background2(1:yyy-1-counter,:),background(yyy+1:end,:));
                counter=counter+1;
            end
        end
        rand_trial=sort(randi(size(background2,1),1,21));
        background2=background2(rand_trial,:);
        vep_bkgd=background2;
        backgroundM=nanmean(background2,1)./norm_vep;
        [ttf_BKGD,f]=pwelch(backgroundM,1500,[],[],Fs);
        ttf_BKGD=ttf_BKGD';
        f=f';
        Bootstat=bootstrp(100,@nanmean,background,1);
        for yy=1:size(Bootstat,1)
            ttf_bkgd_boot(yy,:)=pwelch(Bootstat(yy,:),1500,[],[],Fs);
        end
        
        ttf_bkgd_boot=sort(ttf_bkgd_boot,1);
        ttf_bkgdCI=ttf_bkgd_boot([5 95],:);
        
        for bb=1:length(A)
            temp=abs(f-A(bb));
            temp2=find(temp==min(temp));
            ttf_bkgd_Fr(bb,:)=sum(ttf_BKGD(:,temp2-1:temp2+1),2);
            ttf_bkgdCI_Fr(bb,:)=sum(ttf_bkgdCI(:,temp2-1:temp2+1),2);
        end
        
        ttf_bkgd_Fr(1,:)=ttf_BKGD(:,temp(1));
        ttf_bkgdCI_Fr(1,:)=ttf_bkgdCI(:,temp(1));
        
        
end