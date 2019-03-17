function [processedVEPdata]=preprocessVEP(parsedVEPdata, varargin)

%% Parse input
p = inputParser;
p.addParameter('TemporalFrequency',[1.625 3.25 7.5 15 30],@isnumeric);
p.addParameter('dur_in_sec',1.5,@isnumeric); % length (in seconds) of VEP epoch to take from each trial
p.addParameter('Fs',2000,@isnumeric);
p.addParameter('normalize1',0,@islogical);
p.addParameter('normalize2',0,@islogical);
p.addParameter('normalize3',0,@islogical);
p.addParameter('plot_all',0,@islogical);
p.parse(varargin{:});

dur_in_freq=p.Results.dur_in_sec*p.Results.Fs;
XX=(1:length(parsedVEPdata(1).vep_Fr))/p.Results.Fs;
%% Remove poor data quality VEP

% temporal frequency trials
for XXX=1:size(parsedVEPdata,2)
    vep_FR=parsedVEPdata(XXX).vep_Fr;
    vep_BKGD=parsedVEPdata(XXX).vep_bkgd;
    for w=1:size(vep_FR,1)
        for x=1:size(vep_FR,2)
            if p.Results.plot_all==1
                figure(15)
                plot(XX,squeeze(vep_FR(w,x,:)))
                title(['Frequency: ' num2str(p.Results.TemporalFrequency(w)) ' , Trial: ' num2str(x)])
                ax=gca;
                ax.YLim=[-0.5 0.5];
                ax.XLim=[0 dur_in_freq/p.Results.Fs];
                ax.TickDir='out';
                ax.Box='off';
                hold off
                pause
            end

            if max(abs(squeeze(vep_FR(w,x,:))))>=0.5 || max(abs(squeeze(vep_FR(w,x,:))))<0.02
                if p.Results.plot_all==1
                   xlabel('bad')
                   pause
                end
                vep_Fr(XXX,w,x,:)=vep_FR(w,x,:).*NaN;
            else
                vep_Fr(XXX,w,x,:)=vep_FR(w,x,:);
            end

        end
    end

    % gray background
    for w=1:size(vep_BKGD,1)
    %     figure(15)
    %     plot(XX,vep_bkgd(w,:))
    %     title(['Background, trial: ' num2str(w)])
    %     ax=gca;
    %     ax.YLim=[-0.5 0.5];
    %     ax.XLim=[0 dur_in_freq/p.Results.Fs];
    %     ax.TickDir='out';
    %     ax.Box='off';
    %     hold off

        if max(abs(vep_BKGD(w,:)))>=0.5 || max(abs(vep_BKGD(w,:)))<0.02
            vep_bkgd(XXX,w,:)=vep_BKGD(w,:).*NaN;
        else
            vep_bkgd(XXX,w,:)=vep_BKGD(w,:);
        end
        end
end

%% Normalize signal - mean to 0
if p.Results.normalize1==1
    epoch=(0:0.5:p.Results.dur_in_sec)*p.Results.Fs;
    temp_all_epoch=[];
    for XXX=size(vep_Fr,1)
        for w=1:size(vep_Fr,2)
            for x=1:size(vep_Fr,3)
                for y=1:length(epoch)-1
                    if y==length(epoch)-1
                        temp_epoch=vep_Fr(XXX,w,x,epoch(y)+1:epoch(y+1)+1);
                    else
                        temp_epoch=vep_Fr(XXX,w,x,epoch(y)+1:epoch(y+1));
                    end
                    temp_mean=nanmean(squeeze(temp_epoch));
                    temp_epoch=temp_epoch-(temp_mean*ones(size(temp_epoch,1),size(temp_epoch,2)));
                    temp_all_epoch=cat(4,temp_all_epoch,temp_epoch);
                end
                if p.Results.plot_all==1
                    plot(XX,squeeze(squeeze(vep_Fr(XXX,w,x,:))),'Color',[0.5 0.5 0.5])
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

                vep_FrN(XXX,w,x,:)=temp_all_epoch;
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
                plot(XX,squeeze(vep_bkgd(XXX,w,:)),'Color',[0.5 0.5 0.5])
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

            vep_bkgdN(XXX,w,x,:)=temp_all_epoch;
            temp_all_epoch=[];
        end

        clear vep_bkgd vep_Fr
        vep_bkgd=vep_bkgdN;
        vep_Fr=vep_FrN;
        clear vep_bkgdN vep_FrN
    end
end

%% Normalize timeseries by subtracting the average across all frequencies
if p.Results.normalize2==1
    
    for XXX=1:size(vep_Fr,1)
        norm_vep=nanmean(nanmean(vep_Fr(XXX,:,:,:),3),2);
        for x=1:size(vep_Fr,2)
            for y=1:size(vep_Fr,3)
                vep_FrN(XXX,x,y,:)=vep_Fr(XXX,x,y,:)-norm_vep;

            if p.Results.plot_all==1
                plot(XX,squeeze(squeeze(vep_Fr(XXX,x,y,:))),'Color',[0.5 0.5 0.5])
                hold on
                plot(XX,squeeze(squeeze(vep_FrN(XXX,x,y,:))),'k')
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

        norm_vep=nanmean(vep_bkgd(XXX,:,:),1);
        for x=1:size(vep_Fr,2)
            vep_bkgdN(XXX,x,:)=vep_bkgd(x,:)-norm_vep;
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
    
end

%% Normalize setting max of the mean 
if p.Results.normalize3==1
        norm_vep=max(max(max(abs(nanmean(vep_Fr,3)))));
        processedVEPdata(1).norm_vep=norm_vep;
        processedVEPdata(2).norm_vep=norm_vep;
        processedVEPdata(3).norm_vep=norm_vep;
end

processedVEPdata(1).vep_Fr=squeeze(vep_Fr(1,:,:,:));
processedVEPdata(2).vep_Fr=squeeze(vep_Fr(2,:,:,:));
processedVEPdata(3).vep_Fr=squeeze(vep_Fr(3,:,:,:));

processedVEPdata(1).vep_bkgd=squeeze(vep_bkgd(1,:,:,:));
processedVEPdata(2).vep_bkgd=squeeze(vep_bkgd(2,:,:,:));
processedVEPdata(3).vep_bkgd=squeeze(vep_bkgd(3,:,:,:));
end