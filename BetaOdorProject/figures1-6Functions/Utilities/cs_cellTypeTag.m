

%Determines firing rate for each cell and tags cellinfo with either 'pyr'
%or 'int' accordingly (cutoff is 8.5 Hz). For edge cases, 8 < FR < 9, uses
%spikewidth, where cutoff is 0.35ms

animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
%topDir = 'D:\OdorPlaceAssociation\';
%topDir = 'F:\Data\OdorPlaceAssociation\';
[topDir] = cs_setPaths();

regions = {'CA1','PFC'};


for r = 1:length(regions)
    region = regions{r};
    
    allFR = [];
    for a = 1:length(animals)
        animal = animals {a};
        
        dataDir = [topDir,animal,'Expt\',animal,'_direct\'];
        cd(dataDir)
        cellinfofile = dir([animal,'cellinfo.mat']);
        load(cellinfofile.name);
        
         
        cellfilter = ['(isequal($area,''',region,'''))'];
        cells = evaluatefilter(cellinfo,cellfilter);
        
        days = unique(cells(:,1));
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            load([dataDir, animal, 'spikes', daystr,'.mat']);
            
            daycells = cells(cells(:,1) == day,:);
            epochs = unique(daycells(:,2));
            noeps = daycells(:,[1,3,4]);
            daycells = unique(noeps, 'rows'); 
            
            for c = 1:size(daycells,1)
                %find average over all epochs
                epochFR = [];
                for ep = 1:length(epochs)
                    epoch = epochs(ep);
                    epochspikes = spikes{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)};
                    if ~isempty(epochspikes)
                        if isfield(epochspikes, 'meanrate')
                        epochFR = [epochFR; epochspikes.meanrate];
                        end
                    end                 
                end
                FR = mean(epochFR);
                allFR = [allFR; FR];
                
               
                for ep = 1:length(epochs)
                    epoch = epochs(ep);
                    %add tag
                        if      (length(cellinfo{daycells(c,1)}{epoch}) >= daycells(c,2)) && ...
                            (length(cellinfo{daycells(c,1)}{epoch}{daycells(c,2)}) >= daycells(c,3)) && ... 
                            ~isempty(cellinfo{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)})
                        if FR <= 8.5
                        cellinfo{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)}.type = 'pyr';
                        else
                        cellinfo{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)}.type = 'int';
                        end
                        
                        %use spike width if too close to cutoff
                        if FR < 9 && FR > 8
                            % JHB catch claire did not add the spikewidth
                            % in a consistent place, there are three ways
                            % to get it:
                            spikewidth=nan;
                            if isfield(spikes{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)},'spikewidth')
                               spikewidth = spikes{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)}.spikewidth;
                            end
                            if isnan(spikewidth) && isfield(cellinfo{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)},'spikewidth')
                               spikewidth = cellinfo{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)}.spikewidth;
                            end
                            if isnan(spikewidth)
                                try
                                % this requires raw data...
                                spikewidth = lhcs_getSpikeWidth(topDir, animal, day, daycells(c,2), daycells(c,3));
                                end
                            end

                            if spikewidth > 0.35 % if its a wide width its a pyram, if not or unknown, int
                                cellinfo{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)}.type = 'pyr';
                            else
                                cellinfo{daycells(c,1)}{epoch}{daycells(c,2)}{daycells(c,3)}.type = 'int'; 
                            end
                        end
                    end
                end
            end
        
        end
    save([topDir,animal,'Expt\',animal,'_direct\',animal,'cellinfo.mat'],'cellinfo');
    
    end
    highFR = find(allFR >= 10);
    allFR(highFR) = 10;
    
    %plot distribution
%     figure,
%     histogram(allFR,[0 1 2 3 4 5 6 7 8 9 10])
%     title(region)
end