
function varargout = process_compute_connectivity_sources( varargin )

% @=============================================================================
% Showed PLV, WPLI and OCC to a specif brain lobe.
% INPUT:
% - Connectivity Method: PLV, WPLI or OCC.
% - Type of coupling: Global, Frontal Lobe, Central Lobe, Parietal Lobe,
%   Occipital Lobe, Temporal Lobe.
% - Number of trials for bootstrap.
% =============================================================================@
%
varargout =  {};
%profile on;
eval(macro_method);
%profile report;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Connectivity Sources(64 channels)';
sProcess.FileTag     = '__';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Connectivity';
sProcess.Index       = 901;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'timefreq'};
sProcess.OutputTypes = {'timefreq'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

sProcess.options.method.Comment = 'Connectivity Method:';
sProcess.options.method.Type    = 'combobox';
sProcess.options.method.Value   = {1, {'plv', 'wpli', 'oCC'}};

sProcess.options.ncat.Comment = 'Number of categories:' ;
sProcess.options.ncat.Type    = 'value' ;
sProcess.options.ncat.Value   = {2,'',0};

sProcess.options.trials.Comment = 'Number of trials:' ;
sProcess.options.trials.Type    = 'value' ;
sProcess.options.trials.Value   = {10,'',0};

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
pBar = bst_progress('start','Estimating Electrodes Coupling','Estimating Electrodes Coupling',0,numel(sInputs));

id_method = sProcess.options.method.Value{1}; % acquisisco il metodo
method = sProcess.options.method.Value{2}{id_method};

ncat = sProcess.options.ncat.Value{1}; % acquisisco il numero di trials per Permutation Test e Bootstrap

nRiffle = sProcess.options.trials.Value{1}; % acquisisco il numero di trials per Permutation Test e Bootstrap

dataNames = {'_obs', '_sur'};
dataSession = {'_1', '_2'};
subjectNames = {sInputs.SubjectName};
subjectNames = unique(subjectNames);
nSubjects = numel(subjectNames);

tmpcell = cell([nSubjects,length(dataNames), length(dataSession)]);

TimeMin = zeros(nSubjects,1);
TimeMax = zeros(nSubjects,1);
for rr = 1:ncat
    for subjectIdx = 1:nSubjects
        
        for dataTypeIdx = 1:numel(dataNames)
             
            subjectMask = ~cellfun(@isempty,regexp({sInputs.FileName},unique(subjectNames(subjectIdx)),'match'));
            typeMask= ~cellfun(@isempty,regexp({sInputs.FileName},dataNames(dataTypeIdx),'match'));
            sessions = ~cellfun(@isempty,regexp({sInputs.Comment},dataSession(rr),'match'));
            fileIdx = find(subjectMask & typeMask & sessions);
            % read one file at time
            data = in_bst_timefreq(sInputs(1,fileIdx).FileName);

            TimeMin(subjectIdx) = data.Time(1);
            TimeMax(subjectIdx) = data.Time(2);
            
            % find bad channel
% %             parent = bst_process('GetInputStruct', data.DataFile);
% %             parentStruct{subjectIdx,:} = parent;
% %             parentData = in_bst(parent.FileName);
   
            tmpcell{subjectIdx,dataTypeIdx, rr} = data.TF; % nSubj x nType x nSessions
        end
    end
end


%% Plot WPLI/OCC/iPLV

obs = squeeze(tmpcell(:,1,:)); % nSubjects x nSessions
sur = squeeze(tmpcell(:,2,:)); % nSubjects x nSessions

compute_bootstrap(obs, sur, method, ncat, nRiffle); % listobs [nSubjx1] --> [nChans^2 x nFreqs]

OutputFiles = [];
end


function [meanCatObs, percCatObs, meanCatSur, percCatSur] = compute_bootstrap(listobs, listsur, method, ncat, ntrials)

nFreqs = size(listobs{1},2);

%% SETUP FOR BOOTSTRAPPING
percs       = [2.5, 97.5];

colors = colormap_ripples(true); % setto la colormap
step = floor(size(colors,1)/(ncat-1));
id_color = 1:(step-1):size(colors,1);
colorsObs = colors(id_color,:);

scrsz = get(0,'ScreenSize'); scrsz(2) = scrsz(3)/2; scrsz = scrsz.*[1 1 1/2 1/2];
figure('Position',scrsz);

meanCatObs  = zeros(ncat,nFreqs);
meanCatSur  = zeros(ncat,nFreqs);

percCatObs  = zeros(length(percs), nFreqs, ncat);
percCatSur  = zeros(length(percs), nFreqs, ncat);


%% BOOTSTRAPPING

for c = 1:ncat
    
    
    listobs_tmp(:,1) = listobs(:,c); % prendo una sessione alla volta per tutti i soggetti
    listsur_tmp(:,1) = listsur(:,c);

    bootMeanObs = bootstrap_subjects_eachlobes(listobs_tmp, ntrials);
    bootMeanSur = bootstrap_subjects_eachlobes(listsur_tmp, ntrials);
    
    meanObs = zeros(size(listobs_tmp,1), nFreqs); 
    meanSur = zeros(size(listsur_tmp,1), nFreqs); 
    for ss = 1:size(listobs_tmp,1) 
        meanObs(ss,:) = mean(listobs_tmp{ss,1},1);
        h{:,c} = plot(meanObs(ss,:), 'Color', colorsObs(c,:), 'Linestyle', ':', 'LineWidth', 0.5); hold on; 
        meanSur(ss,:) = mean(listsur_tmp{ss,1},1);
    end
    meanCatObs(c,:) = mean(meanObs,1);
    meanCatSur(c,:) = mean(meanSur,1);
    
    percCatObs (:,:,c) = prctile(bootMeanObs, percs);
    percCatSur (:,:,c) = prctile(bootMeanSur, percs);
    %         meanCatObs(c,:)   = mean(cell2mat(bootMeanObs(c,:)),2).';
    %         percCatObs(:,:,c) = prctile(cell2mat(bootMeanObs(c,:)).', percs);
    %
    %         meanCatSur(c,:)   = mean(cell2mat(bootMeanSur(c,:)),2).';
    %         percCatSur(:,:,c) = prctile(cell2mat(bootMeanSur(c,:)).', percs);
    % %
    h1{:,c} = plotfillbetween_eachlobes(percCatObs(:,:,c), meanCatObs(c,:), colorsObs(c,:), '-',2); hold on;
    %  plotfillbetween_eachlobes(percCatSur(:,:,c), meanCatSur(c,:), colorsObs{c}, ':'); hold on;
    h2{:,c} = plot(meanCatSur(c,:), 'Color', colorsObs(c,:), 'LineStyle', '-.');
   
end

title(upper(method))

xlabel('Freqs [Hz]')
ylabel(upper(method))
xlab = 1:3:40;
xticklabels(xlab);
legend ([h1{:,1}, h2{:,1}, h{:,1}, h1{:,2}, h2{:,2}, h{:,2}], 'Run1 obs', 'Run1 sur', 'Run1 (subj)', 'Run2 obs', 'Run2 sur', 'Run2 (subj)');
% legend(['Run1 obs'],['Run1 sur'], ['Run2 obs'], ['Run2 sur']);
hold off;
end

function bootcomp = bootstrap_subjects_eachlobes(listsubj,ntrials)

nSubjs = size(listsubj,1);
nFreqs = size(listsubj{1},2);
meanchan = zeros(2,nFreqs);

for i = 1:ntrials
    idxsubj = randi([1 length(listsubj)],1,length(listsubj));
    for j = 1:nSubjs
        listsubj_new{j,1} = listsubj{idxsubj(j)};
        meanchan(j,:) = mean(listsubj_new{j,1}); % medio sulle sources
    end
    bootcomp(i,:) = mean(meanchan);
end

end

function [ph] = plotfillbetween_eachlobes(percCat, meanCat, color, linest, linewidth)
% PLOTFILLBETWEEN Plots bootstrapped mean for two lists of patients
%
%  Obs is a nChans^2 x nFreqs x nSubjs
%
% axes(ax);
nFreqs = size(meanCat,2);

x = 1:nFreqs;
y1 = percCat(1,:,1);
y2 = percCat(2,:,1);

ph = plot(x, meanCat(1,:), 'LineStyle', linest, 'LineWidth', linewidth);
ph.Color = color;

hold on;

fh = fill([x fliplr(x)], [y1 fliplr(y2)], color, 'HandleVisibility', 'off');

fh.EdgeColor = 'none';
fh.FaceAlpha = 0.15;

hold off;
end