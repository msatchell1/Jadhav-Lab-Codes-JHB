
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
regions = {'CA1','PFC'};

[topDir] = cs_setPaths();

dataDir = [topDir,'AnalysesAcrossAnimals\'];


for r = 1:length(regions)
region = regions{r};
    
    interneurons =[];
    for a = 1:length(animals)
        animal = animals{a};
        animDir = [topDir, animal, 'Expt\',animal,'_direct\'];
        
        load([animDir,animal,'cellinfo.mat'])
        
        cellfilter = ['isequal($area,''',region,''') & (strcmp($type,''int'')) & ($numspikes > 100)']; 
        
        animcells = evaluatefilter(cellinfo,cellfilter);
        
        noeps = animcells(:,[1 3 4]);
        animcells = unique(noeps,'rows');
        
        days = unique(animcells(:,1));
       
        animvect = repmat(a, size(animcells,1), 1);
        animcells = [animvect, animcells];
        
        interneurons = [interneurons; animcells];
    end
    
    save([dataDir,'interneurons_',region,'.mat'], 'interneurons')
end