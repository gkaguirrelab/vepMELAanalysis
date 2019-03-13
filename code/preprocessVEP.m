function [processedVEPdata]=preprocessVEP(vep_Fr, vep_bkgd, varargin)

%% Parse input
p = inputParser;
p.addParameter('TemporalFrequency',[1.625 3.25 7.5 15 30],@isnumeric);
p.addParameter('dur_in_sec',1.5,@isnumeric); % length (in seconds) of VEP epoch to take from each trial
p.addParameter('Fs',2000,@isnumeric);
p.addParameter('normalize1',0,@islogical);
p.addParameter('normalize2',0,@islogical);
p.addParameter('plot_all',0,@islogical);
p.parse(varargin{:});

dur_in_freq=p.Results.dur_in_sec*p.Results.Fs;
XX=(1:length(vep_Fr))/p.Results.Fs;
%% Remove poor data quality VEP

% temporal frequency trials
for w=1:size(vep_Fr,1)
    for x=1:size(vep_Fr,2)
        if p.Results.plot_all==1
            figure(15)
            plot(XX,squeeze(vep_Fr(w,x,:)))
            title(['Frequency: ' num2str(p.Results.TemporalFrequency(w)) ' , Trial: ' num2str(x)])
            ax=gca;
            ax.YLim=[-0.5 0.5];
            ax.XLim=[0 dur_in_freq/p.Results.Fs];
            ax.TickDir='out';
            ax.Box='off';
            hold off
            pause
        end
        
        if max(abs(squeeze(vep_Fr(w,x,:))))>=0.5 || max(abs(squeeze(vep_Fr(w,x,:))))<0.02
            if p.Results.plot_all==1
               xlabel('bad')
               pause
            end
            vep_Fr(w,x,:)=vep_Fr(w,x,:).*NaN;
        end
        
    end
end

% gray background
for w=1:size(vep_bkgd,1)
%     figure(15)
%     plot(XX,vep_bkgd(w,:))
%     title(['Background, trial: ' num2str(w)])
%     ax=gca;
%     ax.YLim=[-0.5 0.5];
%     ax.XLim=[0 dur_in_freq/p.Results.Fs];
%     ax.TickDir='out';
%     ax.Box='off';
%     hold off

    if max(abs(vep_bkgd(w,:)))>=0.5 || max(abs(vep_bkgd(w,:)))<0.02
        vep_bkgd(w,:)=vep_bkgd(w,:).*NaN;
    end
end


% % Calculate fourier transform w is temporal frequency of the stimuli, x is
% % the repeats (concatonated across sessions)
% for w=1:size(vep_Fr,1)
%     for x=1:size(vep_Fr,2)
%        
%         % Fourier transform
%         ft=fft(squeeze(vep_Fr(w,x,:)));
%         P = abs(ft/p.Results.dur_in_freq);
%         P_data(w,x,:) = P(1:p.Results.dur_in_freq/2+1);
%         
%         
%         % select power spectra for each frequency
%         temp=find(f>=p.Results.TemporalFrequency(w));
%         P_dataFr(w,x)=max(squeeze(P_data(w,x,temp(1,1)-1:temp(1,1))))';
%         clear temp
%     end
%     
% end


%% Normalize signal - mean to 0
if p.Results.normalize1==1
    epoch=(0:0.5:p.Results.dur_in_sec)*p.Results.Fs;
    temp_all_epoch=[];

    for w=1:size(vep_Fr,1)
        for x=1:size(vep_Fr,2)
            for y=1:length(epoch)-1
                if y==length(epoch)-1
                    temp_epoch=vep_Fr(w,x,epoch(y)+1:epoch(y+1)+1);
                else
                    temp_epoch=vep_Fr(w,x,epoch(y)+1:epoch(y+1));
                end
                temp_mean=nanmean(squeeze(temp_epoch));
                temp_epoch=temp_epoch-(temp_mean*ones(size(temp_epoch,1),size(temp_epoch,2)));
                temp_all_epoch=cat(3,temp_all_epoch,temp_epoch);
            end
            if p.Results.plot_all==1
                plot(XX,squeeze(vep_Fr(w,x,:)),'Color',[0.5 0.5 0.5])
                hold on
                plot(XX,squeeze(temp_all_epoch),'k')
                ax=gca;
                ax.YLim=[-0.5 0.5];
                ax.XLim=[0 dur_in_freq/p.Results.Fs];
                ax.TickDir='out';
                ax.Box='off';
                hold off
                pause
            end
            
            vep_FrN(w,x,:)=temp_all_epoch;
            temp_all_epoch=[];
        end
    end

    clear temp_epoch temp_mean
    temp_all_epoch=[];

    for w=1:size(vep_bkgd,1)
        for y=1:length(epoch)-1
            if y==length(epoch)-1
                temp_epoch=vep_bkgd(w,epoch(y)+1:epoch(y+1)+1);
            else
                temp_epoch=vep_bkgd(w,epoch(y)+1:epoch(y+1));
            end
            temp_mean=nanmean(squeeze(temp_epoch));
            temp_epoch=temp_epoch-(temp_mean*ones(size(temp_epoch,1),size(temp_epoch,2)));
            temp_all_epoch=cat(2,temp_all_epoch,temp_epoch);
        end
        if p.Results.plot_all==1
            plot(XX,vep_bkgd(w,:),'Color',[0.5 0.5 0.5])
            hold on
            plot(XX,temp_all_epoch,'k')
            ax=gca;
            ax.YLim=[-0.5 0.5];
            ax.XLim=[0 dur_in_freq/p.Results.Fs];
            ax.TickDir='out';
            ax.Box='off';
            hold off
            pause
        end
        
        vep_bkgdN(w,x,:)=temp_all_epoch;
        temp_all_epoch=[];
    end

    clear vep_bkgd vep_Fr
    vep_bkgd=vep_bkgdN;
    vep_Fr=vep_FrN;
    clear vep_bkgdN vep_FrN
end

%% Normalize timeseries by subtracting the average across all frequencies
if p.Results.normalize2==1
    norm_vep=nanmean(nanmean(vep_Fr,2),1);
    for x=1:size(vep_Fr,1)
        for y=1:size(vep_Fr,2)
            vep_FrN(x,y,:)=vep_Fr(x,y,:)-norm_vep;
        
        if p.Results.plot_all==1
            plot(XX,squeeze(vep_Fr(x,y,:)),'Color',[0.5 0.5 0.5])
            hold on
            plot(XX,squeeze(vep_FrN(x,y,:)),'k')
            ax=gca;
            ax.YLim=[-0.5 0.5];
            ax.XLim=[0 dur_in_freq/p.Results.Fs];
            ax.TickDir='out';
            ax.Box='off';
            hold off
            pause
        end
        end
    end
    
    clear norm_vep;
    
    norm_vep=nanmean(vep_bkgd,1);
    for x=1:size(vep_Fr,1)
        vep_bkgdN(x,:)=vep_bkgd(x,:)-norm_vep;
        if p.Results.plot_all==1
            plot(XX,vep_bkgd(x,:),'Color',[0.5 0.5 0.5])
            hold on
            plot(XX,vep_bkgdN(x,:),'k')
            ax=gca;
            ax.YLim=[-0.5 0.5];
            ax.XLim=[0 dur_in_freq/p.Results.Fs];
            ax.TickDir='out';
            ax.Box='off';
            hold off
            pause
        end
    end
    
    clear vep_bkgd vep_Fr
    vep_bkgd=vep_bkgdN;
    vep_Fr=vep_FrN;
    clear vep_bkgdN vep_FrN
end


processedVEPdata.vep_Fr=vep_Fr;
processedVEPdata.vep_bkgd=vep_bkgd;
end