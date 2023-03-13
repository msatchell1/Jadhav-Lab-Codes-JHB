%% ClaireWrapper/pipeline

% first get your code repository
% cd('C:\Users\Jadhavlab\Documents\gitRepos\Jadhav-Lab-Codes\BetaOdorProject\figure7Functions');
cd('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figure7Functions');
% addpath(genpath('C:\Users\Jadhavlab\Documents\gitRepos\Jadhav-Lab-Codes\BetaOdorProject\figure7Functions'));
addpath(genpath('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figure7Functions'));

% now build your dataset;
usesaved=1;
if usesaved==1
    try % home computer filepath
%         load('G:\BrandeisDatasets\SymanskiData\SuperRat-SymanskiElife-2022-08-02.mat');
        load('D:\OdorPlaceAssociation\SuperRat-SymanskiElife-2022-08-02.mat');
        % remake this so its the folder with this projects codebase
%         cd('E:\GithubCodeRepositories\Jadhav-Lab-Codes\BetaOdorProject');
        cd('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject');

    catch % work computer filepath
        try
            load('E:\BrandeisDatasets\SymanskiData\SuperRat-SymanskiElife-2022-08-02.mat');
            cd('C:\Users\Jadhavlab\Documents\gitRepos\Jadhav-Lab-Codes\BetaOdorProject');
        catch
            [loadfile,filedir]=uigetfile();
            load(fullfile(filedir,loadfile));
            cd('C:\Users\Jadhavlab\Documents\gitRepos\Jadhav-Lab-Codes\BetaOdorProject');
        end
        
    end
else  % or to rebuild the dataset from scratch
    % take the data from the hundreds of files she gave me
    edit ClaireDataAggregator

    % this builds the trajectories (old version)

    % This is the linearization algorithm with four trajectories (2 in, 2 outbound)
    % algorithm:
    % 1. draw a typical trajectory for each of the 4 run types
    %    will ask you to redraw cs39 over again
    %   this is because this is a partial long track, remember no 90* angles on
    %   your tracks!
    % 2. for each run (space between wells) gather the rats distance along and
    %   distance from each trajectory
    % 3. make sure that the odor port entries match which trajectory we think
    %   the rat took (the start port and end port should match whichever
    %   trajectory the rats position remained closest to.
    % 4. Save those data
    %
    % data that are produced:
    % allLinCoords (this contains all the linearized trajectories in a struct)
    % rows: ts, x,y,origin, destination, epoch#, velocity, out left, out right,
    % in left, in right

    % the 4 traj is better if your runs are all wonky (like for short tracks)
    edit LinearizePositionTmazeNTraj
    % first parse pyrams and interneurons based on rate and burst index
    edit ParsePyrIN

end


% some helpful fx for that:
% - one is getcoord_Tmaze, this lets you draw the trajectories, or at
% least is a skeleton for opening a gui and drawing some mean trajectories
%% make single cell examples here:

% to plot some place fields out use this:
edit ClaireQuick2dPlacePlot % i dont think this went into the paper

% plot linearized place fields
edit PlotLinTrajectories

% to calculate the selectivity scores etc for space
edit CalcLinTrajectories %

%% now object coding

% get object coding
edit Object_selectivity
    


%% this is all figure 7 stuff
% get spatial info for object cells
edit PlotObjPlaceCells

edit SummaryObjPlaceInteractions

edit Claire_Odor_Routedecoder

%% and shantanus pca plot

edit sj_save_PCA_plot.m
%% now to analyze ripples...
% one fun question is whether odor representations are replayed with goal
% location representations, although it may be tough as there are only four
% options, so the basrates are very high
edit ClaireRippleAnalysis

%% clear workspace
clearvars -except SuperRat



%% More notes
%{

The main question is this:
do interneurons have higher beta coherence during correct trials?
% there are far fewer error trials

for eeg data
Just pull the eeg from the tetrode with the most cells in that day and use
that

also use the nonreferenced

* we'll use coheren

notes for splitter:
wenbo did a correlation type scoring across the two linearize trajectories
One thing yu can do is run an xcorr between the two rate maps.  The other
way to do this is to bootstrap the absoluite differenes between the rate
maps. the algorithm would go as follows:"
- gather the full ratemaps for the real data, and at each bin take the
absolute difference between the two maps
- to get signficance, gather each 'run' and randomize which trajectory it is.
Then you can build rate maps for each trajectory, and then run the absolute
difference between the two.  WEnbo ran the diff over sum, but
mathematically, that underrepresents high firing pixels and i dont know if
i want that...
the regression i think we're looking for is how similar the rate maps are,
but i'll need to normalize them by maybe their rates overall?



Notes from 4/10/20
for figure 5- put letters for subpanels
-remove odor 1 vs odor 2 rates for cells that are only responsive and not
selective
-add raster for odor responsiveness for some cells
-more example cells
-in supplement, show specifically cells that are odor selective and show
their route selectivity is +, mnopt selective, or opposite selecting
-for pfc, cause the regression works, pick cells that have matching route
selectivity and odor selectivity
-for route selectivity, nspikes per run is fine
-for pfc place fields need to be called spatial fields
-use a binomial test to compare the probability of responsive
cells/selective cells having pfs to those not responsive or selective
-replace the bar graph with pf numbers as x axis with the regression plot
of rates during run and rates during odor presentation
-use claires glm code to decode routes given odor activity.



notes for cross frequency coupling
1. cfc doesnt occur between beta and rr
2. it does occur for some very high and very low freq bands though,
especially in the ob


notes for cell-cell interactions
1. cells have strong oscillations in their cross correlograms
2. some are tightly locked, others less
3. cell pairs are either locked ot beta or rr, generally its the prefrontal
cell that dictates which rhythm (all their ccgrams are that rhythm
4. it kindof appears that for rr CA1 leads, and for beta, PFC leads...
unclear.

5.  it also looks like for beta cells, their best peak is either within 0.005 sec,
or a whole beta cycle away, like 0.035 msec.


%}
% GLM code
%{
glm code
% ok so Cellmatrix is n events by m cells, and its the number of spikes per
event row.  my window will be 800 msec.  So the goal here will be to use
exactly her decoder and then apply it to multiple spatial bin sizes.


% snippit pulled from claires codebase..
K = 5;
cv = cvpartition(numtrials, 'kfold',K);
%        
for k=1:K
    % training/testing indices for this fold
    trainIdx = cv.training(k);
    testIdx = cv.test(k);
    
    % train GLM model
    
    warning('off','all');
    %disp('Fitting GLM')
    mdl = fitglm(Cellmatrix(trainIdx,:), trialIDs(trainIdx),'Distribution', 'binomial');
    
    % predict regression output
    Y_hat = predict(mdl, Cellmatrix(testIdx,:));
    
pred = round(Y_hat); %glm outputs very small values instead of zeros (binomial) for some reason, round to zero
ids = trialIDs(testIdx);
fract_correct(k) = sum(pred == ids)/length(pred);

end
fract_correct = mean(fract_correct);

%Do shuffling
for s = 1:100
    if mod(s,10) == 0 || s == 1
    %disp(['Doing shuffle- iteration ', num2str(s)])
    end
    
    shuff_trialIDs = trialIDs(randperm(length(trialIDs)));
    mdl = fitglm(Cellmatrix(trainIdx,:), shuff_trialIDs(trainIdx),'Distribution', 'binomial');
    Y_hat_shuff = predict(mdl, Cellmatrix(testIdx,:));
    pred_shuff = round(Y_hat_shuff);
    ids_shuff = shuff_trialIDs(testIdx);
    
    fract_correct_shuff(s) = sum(pred_shuff == ids_shuff)/length(pred_shuff);
end
fract_correct_shuff = mean(fract_correct_shuff);

%}
%% random sketchpad

odorinfo=[SuperRat(ses).trialdata.sniffstart SuperRat(ses).trialdata.leftright10,...
    SuperRat(ses).trialdata.CorrIncorr10];
rundata=SuperRat(ses).LinCoords;
% now lets go trial by trial and see what the odor ID is, its correct
% incorrect, and the route the animal took.

% Correct: l is L, 0 is R
%
figure;
for i=1:5
    subplot(5,1,i);
    trial=i;
    runstart=find(rundata(:,1)>odorinfo(trial,1),1,'first');
    runstop=runstart+find(rundata(runstart:end,8)>55,1,'first');
    plot(rundata(runstart:runstop,2),rundata(runstart:runstop,3));
    title(sprintf('odorID=%d, corrincorr=%d',odorinfo(trial,2),odorinfo(trial,3)));
end

%%  This is the beta coherence dataset, this mirrors figure 1 results initially
% coded by CVS

% to calculate beta we use raw lfp files (theyre large)
betaDir='E:\ClaireData';
edit addTetInfo % get the info from tets (especially for ob tets)

%edit gatherEEGdata % pull the eeg data and parse into bands (and verify, old version)
edit gatherEEG % this is the current stable copy of that function

% these figure panels are nolonger in the ms
edit ClaireLFPspatialPlots
% spike-lfp interactions
edit ClaireBetaAnalysis % get spike-band coherence
% lfp only interactions
edit ClaireCrossCoherence % run cross coherence across bands
% also nolonger in dataset
edit ClairePhaseOffset % followup for reviewers
% cell pair analysis, currently only cross-coherence between pfc-ca1 cells,
% and we run it based off of an old paper
edit ClaireCellCellInteractions
