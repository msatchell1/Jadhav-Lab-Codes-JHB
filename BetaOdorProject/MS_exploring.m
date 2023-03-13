% Script for myself to work through and understand the data structure
% SuperRat, as well as the functions used to produce figures in the
% Claire Symanski and John Bladon project. Adapted from John's
% A_Figure7WrapperPipeline.m script.
%
% My code is mostly at the bottom of the script.
%
% ---- Michael Satchell 1/24/23 --------

% First, download John's github repository BetaOdorProject
% at: https://github.com/bladonjay/Jadhav-Lab-Codes-JHB
%
% Then, a copy of the project data is needed that includes the file:
% SuperRat-SymanskiElife-2022-08-02.mat
%
% A good starting point for understanding the code and how to reproduce
% figures is by working through the script John wrote for doing just that:
% Jadhav-Lab-Codes-JHB\BetaOdorProject\figure7Functions\A_Figure7WrapperPipeline.m
% This script imports the SuperRat data structure which is simply a massive
% .mat file containing all the processed cell and behavioral data. I need
% to fully understand this data structure.
% The script also opens Claire's functions used to produce figures. I will
% copy most lines out of A_Figure7WrapperPipeline.m to do something
% similar.

% Add the fig7 code repository to the matlab path
addpath(genpath('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\figure7Functions'));

% Load the SuperRat structure. This is on the hard drive Shantanu gave me.
% NOTE: If loading the original SuperRat structure "SuperRat-SymanskiElife-2022-08-02.mat"
% you will need to then run the functions below that add necessary fields
% to the struct, namely CalcLinTrajctories.m and Object_selectivity.m.
load("E:\OdorPlaceAssociation\SuperRat_added_fields.mat");

% Move working directory to project folder.
cd('C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject');

% Note: I compared the SuperRat structure from John's computer (now on the
% hard drive) with the SuperRat structure on citadel. They appear to be
% identical.

%% make single cell examples here:

% to plot some place fields out use this:
edit ClaireQuick2dPlacePlot % i dont think this went into the paper

% plot linearized place fields. Note this uses SuperRat.LinCoords, a field
% which is empty for sessions 15 - 33. I do not know why.
edit PlotLinTrajectories

% to calculate the selectivity scores etc for space
edit CalcLinTrajectories % Note: This script creates several SuperRat struct fields
% including the follownig:
%     SuperRat.units.PFexist
%     SuperRat.units.FiresDuringRun
%     SuperRat.units.RunRates
%     SuperRat.units.FieldProps
%     SuperRat.units.TrajScores  
%     SuperRat.units.SplitterScore
%
% PFexist is important in being able to run later scripts. For this reason
% I want to save a version of SuperRat that has all this data already in
% it, that way I don't have to run CalcLinTrajectories every time I import
% SuperRat.
%
% NOTE: These new fields are not added for sessions 15-33. This causes
% problems in Object_selectivity() when trying to convert the cell array
% into a matrix. I believe this is because these sessions used a short
% track. 
%
% NOTE: PFexist, FiresDuringRun, RunRates, and LinPlaceFields each have 4
% parts. PFexist and FireDuringRun are 4-element vectors, where each
% element is a logical 1 or 0 in the format: (out left, out right, in left,
% in right) where out refers to rat going to reward well after sniffing
% odor, and in to the return journey.


%% now object coding

% get object coding
edit Object_selectivity % Running this does also update some fields in SuperRat,
% specifically:
%   SuperRat.units.activeCell
%   SuperRat.units.taskResponsive
%   SuperRat.units.OdorRates
%   SuperRat.units.OdorSelective
%   SuperRat.units.startcurves
%   SuperRat.units.endcurves

% I believe that this script deletes the units.OdorMeans field and moves
% that info into units.OdorSelective.

%% Save SuperRat scruct with added fields

% After running the first section in CalcLinTrajectories.m and the sections
% in Object_selectivity.m "Current preprocessing pipeline" and
% "now run odor selectivy calculations (this runs on all blocks as if they were one)"
% the necessary fields will have been added to the SuperRat struct, as well
% as some fields removed. To save this version of the struct for easier use
% run this section.

% filedir = "E:\OdorPlaceAssociation\SuperRat_added_fields.mat";
% save(filedir,"SuperRat",'-v7.3');


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




%% SuperRat Structure Description

% From here all written by Michael Satchell

% Here is what I have been able to understand about the format and
% information of data contained in SuperRat-SymanskiElife-2022-08-02.mat

% SuperRat is a 1x38 struct because there were 38 seperate recording days, a
% variable number per rat ranging between 3 and 7 days. I believe
% that the same cells are recorded across multiple days for each animal,
% although I should check that this is always true.
%

% Fields:
%
% name - name of rat
%
% daynum - recording day
%
% cellinfo - contains info for each cell (on each tetrode as well I
%   believe) such as mean firing rate (Hz?), number of spikes, cell type,
%   etc.
%
% files - contains the directories John used for importing Claire's data
%   into the SuperRat struct. 
%
% oldspikes - format is identical (I think) to cellinfo, except it contains
%   the additional spike data for each cell. 
%
% units - cell and spiking data reformatted by John for easier
%   accessibility.
%
% tracking - Data from tracking rat movement (I think)
%
% trialdata - holds behavioral data related to each trial. 
%
% ripdata - Sharp-wave ripple data.
%
% EpochData - Holds information on epochs. Recordings were done during
%   sleep and running, so the epoch labels indicate this. 
%
% longtrack - indicates if track used had long stem (1) or short stem (0).
% 
% mazeMap - I believe holds position information for maze locations as well
%   as well as some mean beta measurements. "The mazeMap contains the 2d 
%   locations of the goals and odor port as well as the 2d tracts that were
%   used to generate the 1d linearized place fields."
%   Running PlotLinTrajectories.m shows individual run trajectories in 2D
%   space. It appears the homewell is somewhere around (x,y) = (70, 95),
%   and the left reward site is around (27, 15). Note this left reward
%   location seems to shift to around (35, 23) for some sessions...
%
% LinCoords - linear coordinates containing a variety of fields.
%
% AllLinCoords - Similar to LinCoords but with some different fields. Ask
%   about the Lo/Ro and Li/Ri dist fields. ** Ask **
%
% tetinfo - tetrode information. Note that between different days for the
% same animal the tetrode information such as number of cells recorded
% changes.
%
% (The remaining fields simply hold directories to the data used to create
% the SuperRat structure).

% Further Field Explanations:
%
% cellinfo:
%   I believe this data is formatted similar to how Claire's data was
%   organized. For accessing cell and spiking data it may be easier to use 
%   John's reorganized format under the field 'units.' In this format there
%   are multiple layers (dimensions) to the data for each of the 38
%   recordings. A given entry of cellinfo specifies a single rat and recording day.
%   The first dim within that day has between 1 and 9 entries. I believe this
%   would be the number of recordings done on a given day. The second dim
%   is 1x32, and this is the number of tetrodes used in the drives for
%   these recordings. For each tetrode there is a third dim that has length
%   equal to the number of neurons recorded by that tetrode. Each element
%   in this third dim is a 1x1 struct containing the information for that
%   cell with the following fields:
%         {'spikewidth'}
%         {'meanrate'  }
%         {'numspikes' }
%         {'csi'       }
%         {'propbursts'}
%         {'tag'       }
%         {'area'      }
%         {'type'      }
%
% oldspikes - Format seems identical to cellinfo. 
%   John's Description: "spikes just in the cell in a cell method
%   that filterframework used- I just concatenated all the data that was deep
%   in those cell vectors for each cell to get the 'spikes' or 'ts' field."
%   An extra field 'data' is present that holds cell spiking data. This
%   data is presented as a 7-column vector, and I don't yet know what each
%   column is. 
% 
% units - This field is probably the primary field to be used for
%   accessing the experiment data. Each struct under 'units' is a
%   1 x (num_cells) struct that itself has many fields. Each element of the
%   struct represents a cell and the fields provide information on that
%   cell. Field names and descriptions are as follows:
%     {'tet'            } Tetrode number that cell is recorded from.
%     {'unitnum'        } For each tetrode, an ID number for its units. Not
%     continuous because some cells were thrown out (?). 
%     {'ts'             } The cell spiking data. Contains seven columns.
%     From John: "ignore all columns except for the first." The first
%     column is the spike times for that unit.
%     {'meanrate'       } Mean cell firing rate.
%     {'tag'            } 'accepted' or 'mua'. John does not know what this
%     means.
%     {'area'           } Electrode recording area. 'PFC' or 'CA1'.
%     {'type'           } Cell type. 'pyr' or 'in'.
%     {'startcurves'    } "The odor-locked tuning curves locked to odor 
%     start or odor end.  They are the bits that created figure 2abc and d."
%     {'endcurves'      } As mentioned above.
%     {'betamean'       } All the beta and resp fields below are related to
%     the beta and respiratory oscillations measured with LFP electrodes. I
%     don't think I need this information for my analysis, but I should
%     still ask.
%     {'betaMVL'        }
%     {'betaOstat'      }
%     {'betaRstat'      }
%     {'respmean'       }
%     {'respMVL'        }
%     {'respOstat'      }
%     {'respRstat'      }
%     {'betaPhaseDprime'}
%     {'respPhaseDprime'}
%     {'activeCell'     } Logical 1 or 0. 1 if a cell is active (fires at
%     least 100 spikes over all running epochs) and 0 if not. If 0, all the beta and resp
%     fields have NaNs.
%     {'taskResponsive' } 3-element vector. Index 1 is the firing rate
%     during odor decision period, index 2 is firing rate before odor
%     decision period, and index 3 is a p-value for being task (odor)
%     responsive. If p-val is less than 0.05, the cell is considered responsive.
%     {'OdorRates'      } n x 3 matrix, index 1 is firing rate, index 2 is 1 or
%     0, 1 for left odor and 2 for right odor, and index 3 is also 1 or 0,
%     1 for correct choice and 0 for incorrect choice.
%     {'OdorMeans'      } 2 x 2 matrix, all zeros for cells with activeCell
%     = 0. This field is deleted after running Object_selectivity.m
%     {'OdorSelective'  } 2 x 4 matrix, index 1 is diff in mean rates,
%     2 is pval, 3 is p<05 y/n, 4 is the dprime effect size". From reading the
%     code in Object_selectivity.m, the first field 'score' in
%     OdorSelective is the selectivity index mentioned in the section
%     below and in the paper. There are two rows, top is for correct trials
%     and bottom is for incorrect. After running Object_selectivity.m two
%     more columns are added, the mean firing rates for both odors.
%     {'partners'       } Only exists for activeCell = 1. It's basically a 
%     list of cells on other tetrodes whose x correllogram has a significant peak."
%     --------
%     ADD ONS: The following fields are added after running
%     CalcLinTrajectories.m:
%     {'PFexist'        } 4-elem vector, each elem is logical 1 or 0 determining
%       if cell has a place field for trajectories (left out, right out, left in, right in). 
%     {'FiresDuringRun' }
%     {'RunRates'       }
%     {'LinPlaceFields' }
%     {'FieldProps'     }
%     {'TrajScore'      }
%     {'SplitterScore'  }
%
%
% tracking - I think holds data from tracking rat movement with a camera.
%   Has the following fields:
%     {'data'    } 6 column matrix with all the data
%     {'fields'  } descriptions of data type in each column in 'data'.
%     {'filt'    }
%     {'descript'}
% 
% trialdata - Data split up by trial. All times I think are in seconds.
%   Contains following fields:
%     {'starttime'   } Trial start time
%     {'endtime'     } Trial end time 
%     {'leftright10' } Indicates which arm the rat ran to. 1 is for left, 0
%     is for right.
%     {'CorrIncorr10'} Indicates whether the rats choice of arm was correct
%     (1) or incorrect (0).
%     {'sniffstart'  } Time sniffing begins at odor well.
%     {'sniffend'    } Time sniffing ends at odor well.
%     {'rewardstart' } Start time of receiving reward
%     {'rewardend'   } End time of receiving reward
%     {'stem_inTime' } Start time when rat begins to run along stem after
%     sniffing odor? ** Ask **
%     {'stem_outTim' } Time when rat exits stem and enters left/right arm.
%     {'EpochInds'   } 2 columns, first contains same data at 'starttime',
%     and second contains an integer that must indicate epoch number. My
%     guess is that trials that were run back to back belong to a single
%     epoch.
%
% ripdata - struct with many subfields. A mix of Claire's data and John's
%   editing on that data. Important fields:
%    {'starttime'} Ripple start time.
%    {'endtime'} Ripple end time.
%    {'baseline'} Baseline current for SWR in micro amps.
%    {'peak'} Peak current in micro amps.
%    {'tetrode'} Which tetrode SWR was recorded on.
%
% EpochData - struct with length equal to number of epochs. Each substruct
%   contains a field 'type' which is the only field if 'type' = "sleep". If
%   'type' = "run" then there are two additional fields
%     {'type'       } Animal behavior during epoch. "sleep" or "run".
%     {'linearcoord'} Some sort of coordinates that is likely related to
%     the SuperRat.LinCoords field, but I don't know how.
%     {'environment'} Environment behavior occured in. Seems to only be
%     "odorplace".
%
% LinCoords - n x 15 struct that only exists for long track recordings.
%   Subfields are:
%     {'ts'             } Different from the similar field units.ts. This 
%     field is here for: "The camera is strobing at a different rate than 
%     any of the digital acquisitions- its about 30 hz.  So those are the
%     timestamps of each image frame based off of the 'master clock.'"
%     {'x'              } John: "pixelx". I think this means the pixel
%     location in the x-coordinate for the led light atop the rats head
%     used for tracking.
%     {'y'              } John: "pixely average value across 2 LEDs".
%     {'originWell'     } Indicates where the rat is running from. Indexed
%     from 1 to 3, indicating the odor well and 2 reward locations. The
%     odor well = 3, right reward well = 2, and left reward well = 1.
%     {'destinationWell'} Indicates where rat is running to in same way as
%     'originWell'.
%     {'epoch'          } Epoch number.
%     {'speed'          } Rat speed in cm/s
%     {'prefTrajDist'   } This holds the linearized distance travelled from
%     origin to destination.
%     {'prefTrajInd'    } John: "I believe this is horizontal distance from
%     the most stereotypic path".
%
% 


%% Additional Notes

% Cell type clarification:
% 
% Total number of recorded cells:
% CA1 = 1,309; PFC = 717
%
% Active Cells - fire at least 100 spikes across all running epochs of all
% trials. CA1 = 934 (813 (87%) pyramidal, 121, (13%) interneurons);
% PFC = 508 (464 (91%) pyramidal, 44 (9%) interneurons)
%
% Odor Period Active Cells - have to fire at least as many spikes as there
% are trials across all odor periods. CA1 = 170/813, 21%; PFC = 234/464, 50%
%
% Task Responsive Cells - "exhibit significant changes in firing rate following the onset of odor
% sampling." CA1 = 138/170, 81%; PFC = 185/234, 79%
% 
% Choice Selective Cells - Preferred one behavioral choice over the other,
% regardless of whether that choice is right or wrong.
% "Selectivity was calculated by generating a selectivity index (SI), in which the
% difference between the average firing rate response to each odor on correct trials was divided by the
% sum of the two responses, giving a value between –1 and 1." 
% CA1 = 47/138, 34%; PFC = 59/185, 31%
%
% Splitter Cells - "Hippocampal Place Cells Whose Firing Is Modulated by 
% Where the Animal Is Going or Where It Has Been." In this experiment I
% think this corresponds to place cells along the stem of the track that
% fire at different rates depending on which arm the rat is going to
% choose.



%% The first plotting section of Object_Selectivity.m is giving me an error
% when trying to convert the SuperRat units into a matrix. It's telling me
% the field names are not identical between all units, so I am
% investigating that here.

for i = 1:length(SuperRat)
    length(fieldnames(SuperRat(i).units))
end

% The lengths are not the same. This is because sessions 15-33 used a short
% track instead of a long track, and the short track runs are excluded from
% analysis in the portion of the Object_selectivity() code that adds
% additional fields (like units.PFexist).


%% Figuring out how many cells are odor (choice) selective

% % Choice selective and odor selective are used synonymously in the code and
% % paper. Even though some cells are called odor selective, they are really
% % responding to a choice the rat has made. The following will identify
% % the number of choice selective cells for correct and incorrect choices.
% 
% issig_cells = [];
% 
% for ses = 1:length(SuperRat)
%     for i = 1:length(SuperRat(ses).units)
% 
%         OdorSelective = SuperRat(ses).units(i).OdorSelective;
%         if (OdorSelective{1,3} == 1 || OdorSelective{2,3} == 1) &&... % Is the cell odor responsive
%                 SuperRat(ses).units(i).area == "CA1" &&... % What brain area is the cell in
%                 SuperRat(ses).units(i).type == "pyr" &&... % pyr or in
%                 SuperRat(ses).units(i).taskResponsive(3) < 0.05 &&... % Is the cell task responsive
%                 SuperRat(ses).longTrack == 1 % does the track have a long stem
%                 %any(SuperRat(ses).units(i).PFexist(1:2) == 1)
%             % if all these are true, add the index i to this list to count
%             % the number of odor selective cells.
%             issig_cells(end+1) = i;
%         end
%     end
% end
% 
% length(issig_cells) 
% % I get 71 odor selective PFC neurons, and 61 odor selective CA1 neurons.
% % These are only pyramidal task responsive cells as well, but this is more
% % than what is mentioned in the paper. What can account for this?
% %       - Maybe because I am looking at long and short track? Its possible
% %       the paper only considers long track... when restricting to only
% %       long track, I get 46 PFC pyr cells, and 43 CA1 pyr cells.
% 
% 
% % The following is code that John uses to sort out odor selective cells or
% % cells with place fields. This bit comes from PlotObjPlaceCells.m:
% % if contains(Params.nameappend,'Odor','IgnoreCase',true)
% % CodingCell=cellfun(@(a) a(3)==1, {SuperRat(ses).units.OdorSelective});
% % elseif contains(Params.nameappend,'Place','IgnoreCase',true)
% % CodingCell=cellfun(@(a) any(a(1:2)==1), {SuperRat(ses).units.PFexist}); 
% %
% % And this is how he selects for only pyramidal cells in a certain region:
% % inRegion=cellfun(@(a) contains(a,Params.region), {SuperRat(ses).units.area});
% % isPyram=cellfun(@(a) contains(a,'pyr'),{SuperRat(ses).units.type}); % already classified
% % UseCells=find(CodingCell & inRegion & isPyram);
% 
% 
% % The following code does the same thing, but instead comes from
% % Object_selectivity.m. I think this is the way to do it with the updated
% % SuperRat structure that I am using:
% % regCells=allCells(strcmpi({allCells.area},regions{i}));
% % responseTable.allTot(i)=length(regCells); % all cells
% % activeCells=regCells(cell2mat({regCells.activeCell})==1);
% % responseTable.allActive(i)=length(activeCells);
% % 
% % pyrams=activeCells(strcmpi({activeCells.type},'pyr'));
% % responseTable.pyrTot(i)=length(pyrams); % all pyrams
% % % active that are responsive
% % myResp=cellfun(@(a) a(3)<pcrit, {pyrams.taskResponsive});
% % responseTable.pyrResp(i)=sum(myResp);
% % 
% % % now % of task responsive that are odor selective
% % responseTable.pyrSel(i)=sum(cellfun(@(a) a{1,3}==1, {pyrams(myResp).OdorSelective}));


%% Assign Cell IDs
% % The SuperRat struct doesn't seem to include a field for cell ID under
% % units, which I think will be useful. So, I add an ID here for each cell
% % with the format: 'name-day#-tetrode#-unit#', where name is the rat name,
% % day# is the recording day for that rat, tetrode# is the tetrode this
% % unit is recorded on and listed under the units.'tet' field, and unit# is listed
% % under units.'unitnum.'
% 
% 
% for ses = 1:length(SuperRat)
%     for cell = 1:length(SuperRat(ses).units)
% 
%         unitID = strcat(SuperRat(ses).name, "-" ,num2str(SuperRat(ses).daynum), "-",...
%             num2str(SuperRat(ses).units(cell).tet), "-", num2str(SuperRat(ses).units(cell).unitnum));
% 
%         SuperRat(ses).units(cell).unitID = unitID; % Create new field and assign unit ID
% 
%     end
% end
% 
% % Note this cell ID field is not yet used.

%% Reports the number of cells after applying filters

% Parameters for cells to select. These only determine HOW to sort the
% cells, not WHETHER to sort. If you don't want to sort based on any number
% of these fields, just comment out the corresponding code below.
Params.regions = {'PFC'}; % PFC, CA1, or both.
Params.types = {'pyr'}; % pyr, in, or both.
Params.active = {1}; % 1, 0 or both. 1 to sort for active cells, 0 to sort for inactive.
Params.taskResponsive = {1,0}; % 1, 0 or both.
Params.odorSelective = {1}; % 1, 0 or both.

% Long track sessions:
LT_boo = extractfield(SuperRat,'longTrack');
LT_inds = find(LT_boo);
% LTcells is a combination of all units from all long track sessions.
Cells = cell2mat({SuperRat(LT_inds).units});


% To select cells from n number of regions, define a matrix of
% logical values where each row is the logical vector for cells being in
% a region. Combine multiple rows for multiple regions.
regions_boo_mat = zeros(length(Params.regions),length(Cells));
for i = 1:length(Params.regions)
    regions_boo_mat(i,:) = strcmpi({Cells.area}, Params.regions{i});
end
% Then compare all rows together with logical OR by using the any()
% function along the first dimension
regions_boo = any(regions_boo_mat, 1);
Cells = Cells(regions_boo); % Select cells from regions.

% Similar process for cell types.
types_boo_mat = zeros(length(Params.types),length(Cells));
for i = 1:length(Params.types)
    types_boo_mat(i,:) = strcmpi({Cells.type}, Params.types{i});
end
types_boo = any(types_boo_mat, 1);
Cells = Cells(types_boo); % Select cell types.

% Similar process for active cells.
active_boo_mat = zeros(length(Params.active),length(Cells));
for i = 1:length(Params.active)
    active_boo_mat(i,:) = [Cells.activeCell] == Params.active{i};
end
active_boo = any(active_boo_mat, 1);
Cells = Cells(active_boo); % Select active/inactive cells.

% Similar process for task responsiveness.
TR_boo_mat = zeros(length(Params.taskResponsive),length(Cells));
for i = 1:length(Params.taskResponsive)
    TR_boo_mat(i,:) = cellfun(@(a) a(3)<0.05, {Cells.taskResponsive}) == Params.taskResponsive{i};
end
TR_boo = any(TR_boo_mat, 1);
Cells = Cells(TR_boo); % Select based on task responsiveness

% Selecting for choice (odor) selectivity as done in John's code, where he
% only considers the row 'correct' in the table to determine if a cell is
% choice selective. Done in the same way as the filters above.
OS_boo_mat = zeros(length(Params.odorSelective),length(Cells));
for i = 1:length(Params.odorSelective)
    OS_boo_mat(i,:) = cellfun(@(a) a{1,3} == 1, {Cells.OdorSelective}) == Params.odorSelective{i};
    %OS_boo_mat(i,:) = cellfun(@(a) any(a{:,3} == 1), {Cells.OdorSelective}) == Params.odorSelective{i};
end
OS_boo = any(OS_boo_mat, 1);
Cells = Cells(OS_boo); % Select based on odor selectivity.

% Print results
fprintf("%d cells found using filters: Regions: %s, Types: %s, Active: %s," + ...
    " Task Responsive: %s, Odor Selective: %s. \n", length(Cells), cell2mat(Params.regions),...
    cell2mat(Params.types), mat2str(cell2mat(Params.active)), mat2str(cell2mat(Params.taskResponsive)),...
    mat2str(cell2mat(Params.odorSelective)));


% Sorts out only those cells that have place fields. I can
% change which trajectories I want to consider by altering which elements
% in PFexist I consider. PFexist holds the information for existance of place fields. It is a
% 4-element vector with a 1 or 0 for place fields on trajectories (out
% left, out right, in left, in right).
Cells = Cells(cellfun(@(a) any(a(1:4)==1), {Cells.PFexist}));

fprintf("%d of above cells have place fields for considered trajectories. \n", length(Cells));


% choice selective = odor selective in the terminology used in the paper
% and this code. However I will call the cells I am sorting out below
% choice selective instead of odor selective to be clear that I am
% separating them based on their firing responses to the future choice they
% will make, not the odor they receive. This must be done after sorting for
% odor selective cells. This is done only considering correct trials.
leftCS_boo = cellfun(@(a) a{1,5} > a{1,6}, {Cells.OdorSelective}); % mean firing rate 
% is greater for left choices.
rightCS_boo =  cellfun(@(a) a{1,6} > a{1,5}, {Cells.OdorSelective}); % mean firing rate 
% is greater for right choices.
leftCScells = Cells(leftCS_boo);
rightCScells = Cells(rightCS_boo);
fprintf("%d of above cells are correct left choice selective. \n" + ...
    "%d of above cells are correct right choice selective. \n" + ...
    "------------------------- \n", length(leftCScells), length(rightCScells));




%% Plotting place fields

% Now I am going to take these choice selective cells and plot their
% individual place fields.

% Parameters for cells to select. These only determine HOW to sort the
% cells, not WHETHER to sort. If you don't want to sort based on any number
% of these fields, just comment out the corresponding code below, although
% this should only be necessary for PFTrajs. 
Params.regions = {'PFC'}; % PFC, CA1, or both.
Params.types = {'pyr'}; % pyr, in, or both.
Params.active = {1}; % 1, 0 or both. 1 to sort for active cells, 0 to sort for inactive.
Params.taskResponsive = {1,0}; % 1, 0 or both. 1 for TR, 0 for non-TR.
Params.odorSelective = {1}; % 1, 0 or both. 1 for OS, 0 for non-OS.
Params.PFTrajs = {1,2,3,4}; % Place field trajectories to consider.
% 1 = out left, 2 = out right, 3 = in left, 4 = in right. Can be any
% combination of these indices.
Params.leftrightCS = {1}; % 1, 0 or both. 1 = left choice selective, 0 = 
% right choice selective. Only use when Params.odorSelective = {1}. 

% Parameters for filtering tracking data
Params.timesmooth = 8; % number of bins to smooth with (check out SmoothMat2).
Params.speedthresh = 3; % Speed threshold in cm/s for considering in place field creation.
Params.maxtimejump = 1; % Maximum time in seconds allowed between running instances.
% Spiking during gaps longer than this will not be considered in place field
% creation. 

% Param for plotting individual place fields. 1 to plot individual fields 
% and pause after each one, 0 to skip plotting.
Params.plotIndvPFs = 0;
Params.savefigs = 0; % 1 to save all generated figs, 0 not to.
% Directory to save files and folders to. If directory DNE, it will be
% created. Note this will need to be adjusted to the regions and leftrightCS params above. 
Params.savedir = "C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\MS_Figs\Place Fields\PFC\left_choice";


% To find the average place fields of all neurons across sessions, I need
% vectors to store the run spike data and run tracking data for passing to
% cell_SmoothPlacePlot.m.
all_RSD = []; % All run spike data.
all_RTD = []; % All run tracking data.

% For averaging the final color maps of each cell.
all_FCM = zeros([89,89,3,0]); % All final color maps.
all_ratemap = zeros([89,89,0]); % All ratemaps.


% Note that there may be problems with correct/incorrect trial labelling
% for animals CS31 and CS33. If I wish to exclude these from analysis,
% start from session 8.
for ses = 1:length(SuperRat)
    
    if SuperRat(ses).longTrack % Only consider long track sessions

        % Get boolean arrays for each condition

        regions_boo_mat = zeros(length(Params.regions),length(SuperRat(ses).units));
        for i = 1:length(Params.regions)
            regions_boo_mat(i,:) = strcmpi({SuperRat(ses).units.area}, Params.regions{i});
        end
        regions_boo = any(regions_boo_mat, 1);

        types_boo_mat = zeros(length(Params.types),length(SuperRat(ses).units));
        for i = 1:length(Params.types)
            types_boo_mat(i,:) = strcmpi({SuperRat(ses).units.type}, Params.types{i});
        end
        types_boo = any(types_boo_mat, 1);
        
        % Similar process for active cells.
        active_boo_mat = zeros(length(Params.active),length(SuperRat(ses).units));
        for i = 1:length(Params.active)
            active_boo_mat(i,:) = [SuperRat(ses).units.activeCell] == Params.active{i};
        end
        active_boo = any(active_boo_mat, 1);

        % Similar process for task responsiveness.
        TR_boo_mat = zeros(length(Params.taskResponsive),length(SuperRat(ses).units));
        for i = 1:length(Params.taskResponsive)
            TR_boo_mat(i,:) = cellfun(@(a) a(3)<0.05, {SuperRat(ses).units.taskResponsive}) == Params.taskResponsive{i};
        end
        TR_boo = any(TR_boo_mat, 1);
        
        % Similar process for odor selectivity. Note whether only
        % correct trial OS is considered, or correct and incorrect ( a{1,3}
        % == 1 or any(a{:,3} ==1) ).
        OS_boo_mat = zeros(length(Params.odorSelective),length(SuperRat(ses).units));
        for i = 1:length(Params.odorSelective)
            OS_boo_mat(i,:) = cellfun(@(a) a{1,3} == 1, {SuperRat(ses).units.OdorSelective}) == Params.odorSelective{i};
%             OS_boo_mat(i,:) = cellfun(@(a) any(a{:,3} == 1), {SuperRat(ses).units.OdorSelective}) == Params.odorSelective{i};
        end
        OS_boo = any(OS_boo_mat, 1);
        
        % Similar process but for choosing left and right choice selective.
        % Note this is only for correct choices.
        LR_boo_mat = zeros(length(Params.leftrightCS),length(SuperRat(ses).units));
        for i = 1:length(Params.leftrightCS)
            LR_boo_mat(i,:) = cellfun(@(a) a{1,5} > a{1,6}, {SuperRat(ses).units.OdorSelective}) == Params.leftrightCS{i};
        end
        LR_boo = any(LR_boo_mat, 1);
        
        % For place field existance on a subset of trajectories.
        PF_boo = cellfun(@(a) any(a(cell2mat(Params.PFTrajs))==1), {SuperRat(ses).units.PFexist});

        % Indicies of cells in this session that meet all conditions.
        fltr_cells = find(regions_boo & types_boo & active_boo & TR_boo & OS_boo & LR_boo & PF_boo);


        if ~isempty(fltr_cells) % If any cells meet the criteria

            % Get trial and tracking data
            trialdata = SuperRat(ses).trialdata;
            trackingdata = SuperRat(ses).tracking.data;
            
            % Find where in tracking data the mouse was in a 'run' epoch
            % (as opposed to a 'sleep' epoch).
            runepoch_boo = ismember(trackingdata(:,6), SuperRat(ses).RunEpochs);
            runepoch_trdata = trackingdata(runepoch_boo,:); % Run-only tracking data.

            % Find where in tracking data the mouse was running too slow.
            tooslow_boo = SmoothMat2(runepoch_trdata(:,5),[0 50],Params.timesmooth)<=Params.speedthresh;

            runTrackingData = runepoch_trdata(~tooslow_boo,:); % Removes tracking 
            % data where rat velocity is below threshold.

            % We don't want to include data when the rat is simply standing
            % still, so we use Params.maxtimejump to exclude these
            % periods.
            % For times t in runTrackingData, if a gap in time between t
            % and t+1 is greater than Params.maxtimejump, then the array
            % index of time t is stored in stand_inds.
            stand_inds = find(diff(runTrackingData(:,1))>Params.maxtimejump);

            % Creates a n x 2 array where n is the number of individual
            % running periods. The first column is the start time of
            % running events, and the second column the stop time. 
            % runTrackingData(1,1) is the first run period time, and 
            % runTrackingData(stand_inds+1,1) is the start times of all the
            % following run periods. runTrackingData(stand_inds) are end
            % times of run periods.
            runPeriods = [[runTrackingData(1,1); runTrackingData(stand_inds+1,1)]...
                [runTrackingData(stand_inds); runTrackingData(end,1)]];


            % Calculate rate and color maps
            for i = 1:length(fltr_cells)
                % Grabs spike data from a cell.
                spikedata = SuperRat(ses).units(fltr_cells(i)).ts(:,1);

                disp([SuperRat(ses).units(fltr_cells(i)).OdorSelective{1,5}, SuperRat(ses).units(fltr_cells(i)).OdorSelective{1,6}])

                % Passes the spike data and running period times to EpochCoords
                % to get only spikes that occured during running periods.
                runSpikeData = EpochCoords(spikedata, runPeriods);

                % Saves data to matrices for averaging later.
                all_RSD = [all_RSD; runSpikeData];
                all_RTD = [all_RTD; runTrackingData];
                
                % To be passed to the function (cell_SmoothPlacePlot)
                % that creates the place field, the data needs to be in 
                % structs with these field names.
                runTrackingStruct = struct('edit_coords', runTrackingData);
                runSpikesStruct = struct('ts', runSpikeData);
                
                % Get the place field
                [ratemap,~,finalcolormap]=cell_SmoothPlacePlot(runTrackingStruct,runSpikesStruct,...
                    'Factor',2,'suppress',1,'gaussdev',1.5,'ColorScheme','parula'); 
                
                % Store the color map for later use.
                all_FCM(:,:,:,end+1) = finalcolormap;
                all_ratemap(:,:,end+1) = ratemap;

                if Params.plotIndvPFs % Plot individual place fields.
                    F = figure();
                    image(finalcolormap); set(gca,'YDir','normal');
                    title(sprintf('Ses %d Unit %d | Max Rate: %.1f Hz',ses,fltr_cells(i),max(linearize(ratemap))));
%                     c = colorbar('Ticks',[0,1], 'TickLabels',[0,round(max(linearize(ratemap)),1)]);
%                     c.Label.String = "Firing Rate (Hz)";
                    box off; axis off;

                    if Params.savefigs
                        filename = sprintf("ses%d_unit%d.jpeg", ses, fltr_cells(i));
                        figdir = Params.savedir + "/PF_Indv/"; % Full dir for figure

                        if not(isfolder(figdir)) % Makes dir if DNE
                            mkdir(figdir)
                        end
                        saveas(F,figdir+filename, 'jpeg');
                    end

                    %pause
                end
            
            end
        end
    end
end


% Take the total spiking data from all filtered cells and create a place
% field. Note the method I use below is not great, because when the spike
% or tracking data have timesteps that appear multiple times, as will
% happen when cells are recorded using similar reference clocks, I use the
% unique() function to remove the additional instances. This doesn't remove
% much of the spiking data, but it removes A LOT of the tracking data. I'm
% not sure how big of a problem this is.
%%%
% One possible fix to this could be to just shift the spike/tracking times
% by a tiny bit so that they are not identical to another time in the
% array.
%%%
% I think there is also a bigger problem with this method, that being that
% spikes are probably matched to location using time, so that having
% multiple time series in the data may cause confusion for where spikes
% occur. I don't think this is fixable, and I will probably have to find a
% way to combine the place fields after they have been generated. 
%%% 
% When I compare the averged ratemap imagesc plot to the plot generated by
% feeding the combined spike data to cell_SMoothPlacePlot.m, the plots are
% somewhat similar, but have very noticable differences.

% Eliminating repeated time stamps without sorting. Sorting vs not sorting
% changes the place fields, although I am not sure either is okay...
tot_RSD_unq = unique(all_RSD,'stable');
[~,RTD_unq_inds,~] = unique(all_RTD(:,1),'stable'); % Get inds based solely on times.
tot_RTD_unq = all_RTD(RTD_unq_inds,:);

% Put data into structs for input to cell_SmoothPlacePlot
tot_runTrackingStruct = struct('edit_coords', tot_RTD_unq);
tot_runSpikesStruct = struct('ts', tot_RSD_unq);

% Get the place field
[tot_ratemap,~,tot_finalcolormap]=cell_SmoothPlacePlot(tot_runTrackingStruct,tot_runSpikesStruct,...
    'Factor',2,'suppress',1,'gaussdev',1.5,'ColorScheme','parula'); 

F1 = figure();
image(tot_finalcolormap); set(gca,'YDir','normal');
title(sprintf("Combined Spike Data Across %d Units \n Regions: %s, Types: %s, Active: %s," + ...
    " Task Responsive: %s, Odor Selective: %s, leftrightCS: %s \n PF Trajs: %s, Max Firing Freq: %.2f Hz", size(all_FCM,4), cell2mat(Params.regions),...
    cell2mat(Params.types), mat2str(cell2mat(Params.active)), mat2str(cell2mat(Params.taskResponsive)),...
    mat2str(cell2mat(Params.odorSelective)), mat2str(cell2mat(Params.leftrightCS)),...
    mat2str(cell2mat(Params.PFTrajs)), max(linearize(tot_ratemap))));


% Averaging the color maps from indv units (NOT GOOD)
% avg_FCM = mean(tot_FCM, 4);
% figure();
% imagesc(avg_FCM); set(gca,'YDir','normal');
% title(sprintf("Averaged Color Maps Across %d Units \n Regions: %s, Types: %s, Active: %s," + ...
%     " Task Responsive: %s, Odor Selective: %s, leftrightCS: %s \n PF Trajs: %s, Max Firing Freq: %.2f Hz", size(tot_FCM,4), cell2mat(Params.regions),...
%     cell2mat(Params.types), mat2str(cell2mat(Params.active)), mat2str(cell2mat(Params.taskResponsive)),...
%     mat2str(cell2mat(Params.odorSelective)), mat2str(cell2mat(Params.leftrightCS)),...
%     mat2str(cell2mat(Params.PFTrajs)), max(linearize(tot_ratemap))));



% Plotting average of individual ratemaps.

% I want to include pixels where at least half of the ratemaps have values.
% Loops through 2D-space portion of matrix.
for row = 1:size(all_ratemap,1)
    for col = 1:size(all_ratemap,2)
        
        avg_vec = all_ratemap(row,col,:); % Vector to be averaged later
        numnan = sum(isnan(avg_vec)); % Number of NaN elements
        if numnan > length(avg_vec)/2 % If more than half of elems are NaN
            % Assign all values in that vector to be NaN
            all_ratemap(row,col,:) = NaN([1,length(avg_vec)]);
        end

    end
end

% avg_ratemap = mean(all_ratemap,3,'omitnan'); % Exclude only pixels where all values are NaN.
% F2 = figure();
% imagesc(avg_ratemap); set(gca,'YDir','normal');
% title(sprintf("Averaged Ratemaps Across %d Units \n Regions: %s, Types: %s, Active: %s," + ...
%     " Task Responsive: %s, Odor Selective: %s, leftrightCS: %s \n PF Trajs: %s, Max Firing Freq: %.2f Hz", size(all_ratemap,3), cell2mat(Params.regions),...
%     cell2mat(Params.types), mat2str(cell2mat(Params.active)), mat2str(cell2mat(Params.taskResponsive)),...
%     mat2str(cell2mat(Params.odorSelective)), mat2str(cell2mat(Params.leftrightCS)),...
%     mat2str(cell2mat(Params.PFTrajs)), max(linearize(avg_ratemap))));


% Same thing as the ratemap above, except I normalize the ratemap values to
% between 0 and 1 before combining them. This prevents rate maps with high
% firing rates from overpowering the low firing rate ratemaps.
norm_ratemaps = zeros(size(all_ratemap));
for i = 1:size(all_ratemap, 3) % Loops each ratemap
    norm_ratemaps(:,:,i) = all_ratemap(:,:,i)./max(all_ratemap(:,:,i), [], 'all');
end

avg_norm_ratemap = mean(norm_ratemaps,3,'omitnan'); % Exclude only pixels where all values are NaN.
% For pixels with NaN values (as permitted above), the average is taken
% only considering the non-NaN values.
F3 = figure();
imagesc(avg_norm_ratemap); set(gca,'YDir','normal');
title(sprintf("Averaged Normalized Ratemaps Across %d Units \n Regions: %s, Types: %s, Active: %s," + ...
    " Task Responsive: %s, \n Odor Selective: %s, leftrightCS: %s, PF Trajs: %s", size(all_ratemap,3), cell2mat(Params.regions),...
    cell2mat(Params.types), mat2str(cell2mat(Params.active)), mat2str(cell2mat(Params.taskResponsive)),...
    mat2str(cell2mat(Params.odorSelective)), mat2str(cell2mat(Params.leftrightCS)),...
    mat2str(cell2mat(Params.PFTrajs))));
% box off; axis off;



% Save figures
if Params.savefigs
    filename1 = sprintf("Combined Spike Data Map.jpeg");
%     filename2 = sprintf("Avg Ratemap");
    filename3 = sprintf("Avg Norm Ratemap.jpeg");
    figdir = Params.savedir + "/PF_Avg/"; % Full dir for figure

    if not(isfolder(figdir)) % Makes dir if DNE
        mkdir(figdir)
    end
    saveas(F1, figdir+filename1, 'jpeg');
%     savefig(F2, figdir+filename2);
    saveas(F3, figdir+filename3, 'jpeg');
end



%% Working on the beta phase locking portion - Build filters

% The overall LFP/EEG data processing looks like this:
% 1) The raw LFP data exists in the CSXXExpt folders in the data location.
% I put copies on the hard drive that is used throughout this script. John
% also has the LFP data, and its on citadel.
% 2) Filters need to be created or imported in order to bandpass filter the
% LFP data. Common filters are theta, beta, respiratory rhythm, and ripple
% filters. John's github repo LFP-analysis holds the code needed to make
% this filters, and they are easy to make. filterspecifications.m shows how
% to make them. designeegfilt.m is a function for building the filter
% kernels. filtereeg2 takes the filter and runs it on raw data to produce
% the filtered verison of the LFP data, along with markations denoting the
% phase of the signal that is being analyzed.
% 3) The filtered data is then saved back to the CSXXExpt folders. The
% location of the data is marked in the SuperRat struct with a file path.

% What I want to do is take that filtered, phase-marked LFP data and plot
% the spikes of choice selective neurons on it. The best way to do this is
% probably with a histogram.

% First, I need to be able to filter the raw LFP data. Below I make the
% beta filter and others. The lower and upper bounds on the filter are taken from
% LFP-analysis/JadhavEEGFilter/filterspecifications.m. This script uses 15
% Hz as the low bound and 30 Hz as the high bound for beta.

srate = 1500; % Sample rate of LFP data (Hz).

betakernel = designeegfilt(srate,15,30); % Gets filter kernel. Passes in sample rate and filter lower and upper bounds (Hz).
betafilter.kernel = betakernel; % Saves kernel to struct for input to filtereeg2.m
betafilter.samprate = srate; % sample rate
betafilter.descript = 'Beta filter 15-30 Hz'; % Description. Must be a char vector.
% filterfile is the full directory and filename of the filter to be saved.
filterfile = "C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\betafilter.mat";
save(filterfile, 'betafilter');

respkernel = designeegfilt(srate,7,8);
respfilter.kernel = respkernel;
respfilter.samprate = srate;
respfilter.descript = 'Respiratory filter 7-8 Hz';
filterfile = "C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\respfilter.mat";
save(filterfile, 'respfilter');

slowripplekernel = designeegfilt(srate,100,130);
slowripplefilter.kernel = slowripplekernel;
slowripplefilter.samprate = srate;
slowripplefilter.descript = 'Slow Ripple filter 100-130 Hz';
filterfile = "C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\slowripplefilter.mat";
save(filterfile, 'slowripplefilter');

%% Gather LFP data from Claire's data files

% This is the preprocessing step where the raw LFP data is filtered using
% the filters defined above. The raw LFP data is taken from the CSXXExpt
% folders and run through the filters, returning the filtered data and
% phase markations. This data is then reorganized into 'rawcontinuous',
% 'respcontinuous','betacontinuous','ripplecontinuous' and along with the 
% filters everything is saved in the EEG folder. There is a version saved
% for each region, day and animal. Everything is saved into .mat files
% called: "CS(ANIMALNUM)day(DAYNUM)(REGION)eegdata.mat" and these files are
% saved in theEEG subfolder under each CSXXExpt folder.

% There are a disgusting number of .mat files in each CSXXExpt and I have no
% idea what any of them hold, but the CSXXExpt files I pulled off John's
% computer did NOT already have copies of the eegdata.mat files that
% gatherEEG.m produces, so this functions needs to be run at least once if
% it has never been run on the dataset you are currently using. It takes
% quite a while to run gathrEEG.m on my laptop, around 2 hours. But once
% run once, the data files are added and the location of the data is
% recorded in the SuperRat data structure under the fields '(REGION)eegFile'.
edit gatherEEG.m;

% This function works with the beta LFP data and should be useful for
% learning what to do. Note I need to have the LFP-analysis code to run the
% script. It also adds data (that might already have been there) to
% SuperRat
edit ClaireBetaAnalysis.m;


%% Save SuperRat struct again after updating eegFile directories and 
% running some of ClaireBetaAnalysis.m. I need to dive more into what
% exactly is changed in SuperRat by this script.

% filedir = "E:\OdorPlaceAssociation\SuperRat_added_fields.mat";
% save(filedir,"SuperRat",'-v7.3');

%% Load LFP data from session 1

LFPp.regions = {'CA1'}; % PFC, CA1, or both.


% Here we load in a single LFP file. Note that in the .mat file, there are
% both filters and data structures. In each data structure is the LFP data
% along with other data. The format is such:
% Column 1: reference clock times in seconds
% Column 2: filtered LFP amplitude
% Column 3: instantaneous phase (from a Hilbert transform?). The values
% fall between -pi and pi, where -pi corresponds to the trough of the LFP
% wave.
% Column 4: the envolope used
% Note the raw unfiltered LFP data only has columns 1 and 2.
LFPstruct = matfile(SuperRat(1).([LFPp.regions{1} 'eegFile']));

% Extract LFP data
betaLFPstruct = LFPstruct.betacontinuous;
% rawLFPstruct = LFPstruct.rawcontinuous;
% respLFPstruct = LFPstruct.respcontinuous;
% rippleLFPstruct = LFPstruct.ripplecontinuous;

% We don't need the LFPstruct anymore and its huge so lets clear it
clear 'LFPstruct'

%% Visualize the filtered vs raw LFP data

% Plot portions of the LFP traces by themselves and with phase demarkations
ind_range = 1:1:1000; % Indices to plot

figure();
hold on;
title("Beta Filtered LFP")
plot(betaLFPstruct(ind_range,1), betaLFPstruct(ind_range,2))
plot(betaLFPstruct(ind_range,1), betaLFPstruct(ind_range,3))
plot(betaLFPstruct(ind_range,1), betaLFPstruct(ind_range,4))
legend(["Beta LFP", "Phase clock", "Envelope"])     
xlabel("Time (s)")

% figure();
% hold on;
% title("Raw LFP")
% plot(rawLFPstruct(ind_range,1), rawLFPstruct(ind_range,2))
% legend("Raw LFP")
% xlabel("Time (s)")
% 
% figure();
% hold on;
% title("Respiratory Filtered LFP")
% plot(respLFPstruct(ind_range,1), respLFPstruct(ind_range,2))
% plot(respLFPstruct(ind_range,1), respLFPstruct(ind_range,3))
% legend(["Resp LFP", "Phase clock"])
% xlabel("Time (s)")
% 
% figure();
% hold on;
% title("Ripple Filtered LFP")
% plot(rippleLFPstruct(ind_range,1), rippleLFPstruct(ind_range,2))
% plot(rippleLFPstruct(ind_range,1), rippleLFPstruct(ind_range,3))
% legend(["Ripple LFP", "Phase clock"])
% xlabel("Time (s)")


%% Plot a single neuron's spiking alongside the beta filtered LFP

ind_range = 1:1:5000; % Indices to plot

% Grab a single neuron from the same session and from the same brain
% region as the LFP data
singleunitspikes = SuperRat(1).units(1).ts(:,1);
% min and max times in time window
LFPtimes_minmax = [min(betaLFPstruct(ind_range,1)), max(betaLFPstruct(ind_range,1))];
% Grabs only spikes that occur in this time window
spikes_inwindow = singleunitspikes((LFPtimes_minmax(1) < singleunitspikes) & (singleunitspikes < LFPtimes_minmax(2)));

figure();
hold on;
title("Beta Filtered LFP and Single Neuron Spiketimes")
plot(betaLFPstruct(ind_range,1), betaLFPstruct(ind_range,2))
plot(spikes_inwindow, zeros(length(spikes_inwindow)), 'o')
legend(["Beta LFP", "Spiketimes"])
xlabel("Time (s)")


%% Isolate the decision making period (nose poke in to nose poke out) and
% plot spikes/beta LFP only during that period.

singleunitspikes = SuperRat(1).units(1).ts(:,1);

% This extracts the sniff start and end times for a single session. Note
% that in ClareBetaAnalysis.m, John excludes sniff sessions that are
% shorter than 0.5s or longer than 2.5s. I am not going to do that yet
% unless Shantanu tells me to.
snifftimes=[SuperRat(1).trialdata.sniffstart SuperRat(1).trialdata.sniffend];

% John's code uses a function EpochCoords.m to return the spikes that belong
% to a series of events. I could do that... Or I could try to do it myself
% and compare with EpochCoords.m. Just as a note, EpochCoords can also
% return cellcoords which shows which spike times fall into which event,
% and segments which I believe gives the indices in spiketimes of whenever
% a transition between an event and non-event occurs, in either direction. 
[sniffspikes, segments, cellcoords] = EpochCoords(singleunitspikes, snifftimes); % Using EpochCoords

% Trying to do the same thing myself:
spikeinds = [];
for i = 1:length(snifftimes)
    temp_inds = find((singleunitspikes > snifftimes(i,1)) & (singleunitspikes < snifftimes(i,2)));
    spikeinds = [spikeinds; temp_inds];
end
sniffspikes2 = singleunitspikes(spikeinds); % It works! This is the same as above.


% Anyways, I will probably just use the EpochCoords function in the future.
% I want to plot both the beta LFP and spikes, but only during decision
% periods (sniff times). To get the LFP I think I can just pass it to
% EpochCoords:
sniffLFPtimes = EpochCoords(betaLFPstruct(:,1), snifftimes);
[~,sniffLFPinds,~] = intersect(betaLFPstruct(:,1), sniffLFPtimes); % It seems a bit strange to 
% me the EpochCoords doesn't return the indices of the new coords in terms
% of the oldcoord array, making me have to get the indices this way.

sniffBetaStruct = betaLFPstruct(sniffLFPinds,:); % Sorts out only times during decision period.

figure();
hold on;
title("Decision Periods Only")
plot(sniffLFPtimes, sniffBetaStruct(:,2))
plot(sniffspikes, zeros(length(sniffspikes)), 'o');
legend(["Beta LFP", "Spiketimes"])
xlabel("Time (s)")
        

%% Generate beta phase histograms
% Shantanu wants me to look in PFC and CA1 separately. Look only during
% decision period. I need to do this for all choice selective cells, and
% also for left vs right choice selective to see if any of these groups
% fire preferentially during a specific beta phase.


% Plot, for a single neuron in a single session, a histogram of spikes
% aligned from 0 to 2pi onto the phase of the LFP beta occuring during that
% session. 

% Repeating what has been done above:
singleunitspikes = SuperRat(1).units(1).ts(:,1);
snifftimes=[SuperRat(1).trialdata.sniffstart SuperRat(1).trialdata.sniffend];
sniffspikes = EpochCoords(singleunitspikes, snifftimes);
sniffLFPtimes = EpochCoords(betaLFPstruct(:,1), snifftimes);
[~,sniffLFPinds,~] = intersect(betaLFPstruct(:,1), sniffLFPtimes);
sniffBetaStruct = betaLFPstruct(sniffLFPinds,:);

% Because I need to assign each spike a phase value, I will need to do some
% interpolation because neither the spike nor the LFP data are continuous,
% and they do not line up with each other, so interpolation will be needed
% to fill in the values. What John does in ClaireBetaAnalysis is nearest
% neighbor interpolation, where the value of the nearest neighboring point
% is simply assumed. Examining the data, I think this should be ok because
% the sampling rate is pretty high compared to the frequency of oscillation for
% the beta wave. However I am not certain that when analyzing higher frequency
% oscillations we would be able to get away with this. 

% Just to be consistent, I will do linear interpolation instead of nearest
% neighbor. Below I feed in the LFP times, LFP phase values, and spike
% times. I get out phase values at the spike times.
spikephases = interp1(sniffBetaStruct(:,1), sniffBetaStruct(:,3), sniffspikes, 'linear');

figure();
h = histogram(spikephases, linspace(-pi, pi, 10+ceil(sqrt(length(spikephases)))));
title("Beta Phase Spike Histogram For a Single Neuron")
set(h,'LineStyle','none');  
set(gca,'XTick',[-pi, 0, pi],'XTickLabel',{'-\pi','0','\pi'});
ylabel("Number of Spikes");
xlabel("Beta Phase");

% This is John's code for plotting the histogram
% figure; 
% ha=histogram([spikephases; spikephases+pi*2],linspace(-pi,pi*3,10+ceil(sqrt(length(spikephases)))),...
%     'Normalization','probability');
% set(ha,'LineStyle','none');
% set(gca,'XTick',[pi.*[-1:4]],'XTickLabel',{'-2\pi','','0','','2\pi'});
% ylabel('Probability of Spike'); xlabel('Beta Phase'); box off;



%% MAIN SECTION FOR BETA PHASE
% Repeat the above for all odor selective cells across all sessions

% Remaining LFP params for cells to include in analysis
LFPp.regions = {'PFC'}; % PFC, CA1, or both.
LFPp.types = {'pyr'}; % pyr, in, or both.
LFPp.active = {1}; % 1, 0 or both. 1 to sort for active cells, 0 to sort for inactive.
LFPp.taskResponsive = {1}; % 1, 0 or both. 1 for TR, 0 for non-TR.
LFPp.odorSelective = {1,0}; % 1, 0 or both. 1 for OS, 0 for non-OS.
LFPp.leftrightCS = {1,0}; % 1, 0 or both. 1 = left choice selective, 0 = 
% right choice selective. Only use when Params.odorSelective = {1}. 
LFPp.savefigs = 1;
LFPp.savedir = "C:\Users\micha\OneDrive - brandeis.edu\Brandeis Stuff\Jadhav Lab\Jadhav-Lab-Codes-JHB\BetaOdorProject\MS_Figs";

all_spikephases = {}; % Array to hold phase at which every spike occurs for
% all considered neurons.

% I made a change to only consider correct trials like John does in his
% analysis. This really helps my results, making me think this is necessary
% to do the phase analysis correctly. This would kind of make sense,
% because if incorrect trials have the odor identity swapped, then the
% cells with opposite choice selectivity might be getting included, messing
% up the results.
% I also included only trials between 0.5 and 2.5s in length as John did.
% This did not seem to make much of a difference in the results.

for ses = 1:length(SuperRat)
    
    if SuperRat(ses).longTrack % Only consider long track sessions

        % Get boolean arrays for each condition

        regions_boo_mat = zeros(length(LFPp.regions),length(SuperRat(ses).units));
        for i = 1:length(LFPp.regions)
            regions_boo_mat(i,:) = strcmpi({SuperRat(ses).units.area}, LFPp.regions{i});
        end
        regions_boo = any(regions_boo_mat, 1);

        types_boo_mat = zeros(length(LFPp.types),length(SuperRat(ses).units));
        for i = 1:length(LFPp.types)
            types_boo_mat(i,:) = strcmpi({SuperRat(ses).units.type}, LFPp.types{i});
        end
        types_boo = any(types_boo_mat, 1);
        
        % Similar process for active cells.
        active_boo_mat = zeros(length(LFPp.active),length(SuperRat(ses).units));
        for i = 1:length(LFPp.active)
            active_boo_mat(i,:) = [SuperRat(ses).units.activeCell] == LFPp.active{i};
        end
        active_boo = any(active_boo_mat, 1);

        % Similar process for task responsiveness.
        TR_boo_mat = zeros(length(LFPp.taskResponsive),length(SuperRat(ses).units));
        for i = 1:length(LFPp.taskResponsive)
            TR_boo_mat(i,:) = cellfun(@(a) a(3)<0.05, {SuperRat(ses).units.taskResponsive}) == LFPp.taskResponsive{i};
        end
        TR_boo = any(TR_boo_mat, 1);
        
        % Similar process for odor selectivity. Note whether only
        % correct trial OS is considered, or correct and incorrect ( a{1,3}
        % == 1 or any(a{:,3} ==1) ).
        OS_boo_mat = zeros(length(LFPp.odorSelective),length(SuperRat(ses).units));
        for i = 1:length(LFPp.odorSelective)
            OS_boo_mat(i,:) = cellfun(@(a) a{1,3} == 1, {SuperRat(ses).units.OdorSelective}) == LFPp.odorSelective{i};
%             OS_boo_mat(i,:) = cellfun(@(a) any(a{:,3} == 1), {SuperRat(ses).units.OdorSelective}) == LFPp.odorSelective{i};
        end
        OS_boo = any(OS_boo_mat, 1);
        
        % Similar process but for choosing left and right choice selective.
        % Note this is only for correct choices.
        LR_boo_mat = zeros(length(LFPp.leftrightCS),length(SuperRat(ses).units));
        for i = 1:length(LFPp.leftrightCS)
            LR_boo_mat(i,:) = cellfun(@(a) a{1,5} > a{1,6}, {SuperRat(ses).units.OdorSelective}) == LFPp.leftrightCS{i};
        end
        LR_boo = any(LR_boo_mat, 1);
        

        % Indicies of cells in this session that meet all conditions.
        fltr_cells = find(regions_boo & types_boo & active_boo & TR_boo & OS_boo & LR_boo);


        if ~isempty(fltr_cells) % If any cells meet the criteria
            %%%%%%%%%%
            % Note this only works for ONE region currently. I cannot do
            % multiple regions. 
            %%%%%%%%%%
            % Load in all LFP data. 
            LFPstruct = matfile(SuperRat(ses).([LFPp.regions{1} 'eegFile']));
            % Extract beta LFP data
            betaLFPstruct = LFPstruct.betacontinuous;
            
            % Get decision period start and end times
            snifftimes = [SuperRat(ses).trialdata.sniffstart SuperRat(ses).trialdata.sniffend];
            % We only want correct trials
            trialCorr = logical(SuperRat(ses).trialdata.CorrIncorr10);
            % Only include trials that are correct and between 0.5 and 2.5s
            % long.
            oktrials = trialCorr & diff(snifftimes,1,2)>=.5 & diff(snifftimes,1,2)<2.5;
            % Only sniffs for correct trials
            snifftimesOK = snifftimes(oktrials,:);
            % Get corresponding decision period betaLFP times.
            sniffLFPtimes = EpochCoords(betaLFPstruct(:,1), snifftimesOK);
            % Get indices for decision period
            [~,sniffLFPinds,~] = intersect(betaLFPstruct(:,1), sniffLFPtimes);
            % New structure that now has all beta LFP data for decision
            % periods only.
            sniffBetaStruct = betaLFPstruct(sniffLFPinds,:);


            for i = 1:length(fltr_cells)
                
                % Get spike times of a unit.
                unitspikes = SuperRat(ses).units(i).ts(:,1);
                % Get spikes that occur during decision period.
                sniffspikes = EpochCoords(unitspikes, snifftimesOK);

                if length(sniffspikes) > 30 % Only neurons that spike at least 30 times during decision periods in session
                    % Interpolate to get phase values at each spike time
                    spikephases = interp1(sniffBetaStruct(:,1), sniffBetaStruct(:,3), sniffspikes, 'linear');
    
                    % Append these phases to the total array.
                    all_spikephases{end+1} = spikephases;
    
                    % Plot individual histogram for each unit
                    figh_indv = figure();
%                     h_indv = histogram(spikephases, linspace(-pi, pi, 10+ceil(sqrt(length(spikephases)))));
                     h_indv = histogram([spikephases; spikephases+pi*2],linspace(-pi,pi*3,10+ceil(sqrt(length(spikephases)))));
                    set(h_indv,'LineStyle','none');  
%                     set(gca,'XTick',[-pi, 0, pi],'XTickLabel',{'-\pi','0','\pi'});
                    set(gca,'XTick',[pi.*[-1:4]],'XTickLabel',{'-2\pi','','0','','2\pi'});
                    ylabel("Number of Spikes");
                    xlabel("Beta Phase");
                    title(sprintf("Beta Phase Spike Histogram Ses %d Unit %d \n Regions: %s, Types: %s, Active: %s," + ...
                        " Task Responsive: %s, \n Odor Selective: %s, leftrightCS: %s", ses, fltr_cells(i), cell2mat(LFPp.regions),...
                        cell2mat(LFPp.types), mat2str(cell2mat(LFPp.active)), mat2str(cell2mat(LFPp.taskResponsive)),...
                        mat2str(cell2mat(LFPp.odorSelective)), mat2str(cell2mat(LFPp.leftrightCS))));


                    if LFPp.savefigs
                        filename_indv = sprintf("Beta Spike Histogram Ses %d Unit %d.jpeg", ses, fltr_cells(i));

                        figdir = LFPp.savedir + "/Phase Plots/" + sprintf("Indv Unit Histograms, Regions %s, Types %s, Active %s," + ...
                        " Task Responsive %s, Odor Selective %s, leftrightCS %s", cell2mat(LFPp.regions),...
                        cell2mat(LFPp.types), mat2str(cell2mat(LFPp.active)), mat2str(cell2mat(LFPp.taskResponsive)),...
                        mat2str(cell2mat(LFPp.odorSelective)), mat2str(cell2mat(LFPp.leftrightCS))) + "/"; % Full dir for figure
                    
                        if not(isfolder(figdir)) % Makes dir if DNE
                            mkdir(figdir)
                        end
                        saveas(figh_indv, figdir+filename_indv, 'jpeg');
                    end
    
                end
            end
        end
    end
    fprintf("Session %d completed \n", ses)
end

% Plot cumulative figure.

figh = figure();
% h = histogram(vertcat(all_spikephases{:}), linspace(-pi, pi, 10+ceil(sqrt(length(vertcat(all_spikephases{:}))))));
h = histogram([vertcat(all_spikephases{:}); vertcat(all_spikephases{:})+pi*2],...
    linspace(-pi,pi*3,10+ceil(sqrt(length(vertcat(all_spikephases{:}))))));
set(h,'LineStyle','none', 'FaceColor', 'red');  
% set(gca,'XTick',[-pi, 0, pi],'XTickLabel',{'-\pi','0','\pi'});
set(gca,'XTick',[pi.*[-1:4]],'XTickLabel',{'-2\pi','','0','','2\pi'});
ylabel("Number of Spikes");
xlabel("Beta Phase");
title(sprintf("Beta Phase Spike Histogram Across %d Units \n Regions: %s, Types: %s, Active: %s," + ...
    " Task Responsive: %s, \n Odor Selective: %s, leftrightCS: %s", size(all_spikephases,2), cell2mat(LFPp.regions),...
    cell2mat(LFPp.types), mat2str(cell2mat(LFPp.active)), mat2str(cell2mat(LFPp.taskResponsive)),...
    mat2str(cell2mat(LFPp.odorSelective)), mat2str(cell2mat(LFPp.leftrightCS))));


if LFPp.savefigs
    filename1 = sprintf("Beta Spike Histogram, Regions %s, Types %s, Active %s," + ...
    " Task Responsive %s, Odor Selective %s, leftrightCS %s.jpeg", cell2mat(LFPp.regions),...
    cell2mat(LFPp.types), mat2str(cell2mat(LFPp.active)), mat2str(cell2mat(LFPp.taskResponsive)),...
    mat2str(cell2mat(LFPp.odorSelective)), mat2str(cell2mat(LFPp.leftrightCS)));
    figdir = LFPp.savedir + "/Phase Plots/"; % Full dir for figure

    if not(isfolder(figdir)) % Makes dir if DNE
        mkdir(figdir)
    end
    saveas(figh, figdir+filename1, 'jpeg');
end

