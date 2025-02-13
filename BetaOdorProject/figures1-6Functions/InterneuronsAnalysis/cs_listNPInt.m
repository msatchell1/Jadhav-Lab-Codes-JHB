clear
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};
%regions = {'CA1'};

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];


for r = 1:length(regions)
    region = regions{r};
    
    npInt =[]; activeInt=[];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        % and we arent using the 100 spikes filter because these are
        % interneurons
        cellfilter = ['isequal($area,''',region,''') && isequal($type, ''int'')'];
        %         cellfilter = ['isequal($type, ''pyr'')']
        %         cellfilter = ['isequal($area,''',region,''')']
        
        cells = evaluatefilter(cellinfo,cellfilter);

        % useful filter if you're taking epoch specific data
        %cells=allcells(ismember(allcells(:,[1 2]), runEps,'rows'),:);

        runEps = cs_getRunEpochs(animDir,animal,'odorplace');
        
        days = unique(runEps(:,1));
         noeps = cells(:,[1 3 4]);
        cells = unique(noeps,'rows');
        spikes = loaddatastruct(animDir, animal, 'spikes',days);
        nosepokeWindow = loaddatastruct(animDir, animal, 'nosepokeWindow',days);
        odorTriggers = loaddatastruct(animDir, animal, 'odorTriggers',days);
        
        for d = 1:length(days)
            day = days(d);
            daystr = getTwoDigitNumber(day);
            
            daycells = cells(cells(:,1) == day,:);
            
            
            eps = runEps((runEps(:,1) == day),2);
            
            
            for c = 1:size(daycells,1)
                
                npspikes = [];
                prespikes = [];
                totaltrigs = 0;
                cell = daycells(c,:);
                triallabels = [];
                
                for ep = 1:length(eps)
                    epoch = eps(ep);
                    
                    if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)})
                        if ~isempty(spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data)
                            epspikes = spikes{cell(1)}{epoch}{cell(2)}{cell(3)}.data(:,1);
                            
                            trigs = nosepokeWindow{day}{epoch};
                            totaltrigs = totaltrigs + size(trigs,1);
                            
                            winlength = nosepokeWindow{day}{epoch}(:,2)- nosepokeWindow{day}{epoch}(:,1);
                            pretrigs = [nosepokeWindow{day}{epoch}(:,1)-winlength, nosepokeWindow{day}{epoch}(:,1)];
                            
                            
                            for t = 1:size(trigs,1);
                                
                                %gather spikes
                                winspikes = epspikes(logical(isExcluded(epspikes, trigs(t,:))));
                                npspikes = [npspikes; length(winspikes)];
                                
                                pre = epspikes(logical(isExcluded(epspikes,pretrigs(t,:))));
                                prespikes = [prespikes;length(pre)];
                            end
                            
                            %gather trial labels
                            labels = zeros(size(trigs,1),2);
                            [cr,cl] = cs_getSpecificTrialTypeInds(odorTriggers{day}{epoch});
                            labels([cl;cr],1) = 1;
                            labels(cl,2) = 1;
                            triallabels = [triallabels; labels];
                        end
                        
                    end
                    
                end
                
                if ~isempty(npspikes) || ~isempty(prespikes)
                    correctleft = find(triallabels(:,1) == 1 & triallabels(:,2) == 1);
                    correctright = find(triallabels(:,1) == 1 & triallabels(:,2) == 0);
                    
                    lspikes = [npspikes(correctleft), prespikes(correctleft)];
                    rspikes = [npspikes(correctright), prespikes(correctright)];
                    allspikes = [lspikes; rspikes];
                    p1 = signrank(lspikes(:,1),lspikes(:,2));
                    p2 = signrank(rspikes(:,1),rspikes(:,2));
                    p3 = signrank(allspkes(:,1),allspikes(:,2));
                    
                    if sum(npspikes)/totaltrigs >=1 %at least 1 spike per trial during NP
                       activeInt=[activeInt; a, cell];

                        if p1 < 0.05 || p2 <0.05 || p3 <0.05
                            npInt = [npInt; a, cell];
                        end
                    end
                end
                
               
            end
        end
    end
    
    save([dataDir,'npInt_',region,'.mat'], 'npInt','activeInt')
end
clear