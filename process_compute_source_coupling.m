function varargout = process_compute_source_coupling( varargin )

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
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
sProcess.Comment     = 'Source Coupling';
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
sProcess.options.method.Value   = {1, {'<none>', 'plv', 'wpli', 'oCC'}};

sProcess.options.coupling.Comment = 'Coupling:';
sProcess.options.coupling.Type    = 'combobox';
sProcess.options.coupling.Value   = {1, {'<none>', 'Total', 'Intra-Hemisph.', 'Inter-Hemisph.', 'Intra-PreRol.', 'Intra-PostRol.', 'PreRol-vs-PostRol'}};

sProcess.options.boots.Comment = 'Bootstrap:';
sProcess.options.boots.Type    = 'value';
sProcess.options.boots.Value   = {1000, 'trials', 0};

sProcess.options.perctl.Comment = 'Percentiles for plot:';
sProcess.options.perctl.Type    = 'range';
sProcess.options.perctl.Value   = {[2.50, 97.50], '', 2};

sProcess.options.nameCat.Comment = 'Categories Name: ';
sProcess.options.nameCat.Type    = 'text';
sProcess.options.nameCat.Value   = 'IGE, C, R';

sProcess.options.folder.Comment = 'Destination folder:';
sProcess.options.folder.Type    = 'text';
sProcess.options.folder.Value   = 'ALL_SUBJECTS';

sProcess.options.subfolder.Comment = 'Destination subfolder:';
sProcess.options.subfolder.Type    = 'text';
sProcess.options.subfolder.Value   = 'list_subjects';

% === FREQ BANDS
sProcess.options.label.Comment = '<BR><U><B>Select brain regions:</B></U>:';
sProcess.options.label.Type    = 'label';
sProcess.options.freqbands.Comment = 'Frequency bands for the Hilbert transform:';
sProcess.options.freqbands.Type    = 'groupbands';
%sProcess.options.freqbands.Value   = bst_get('DefaultFreqBands');
sProcess.options.freqbands.Value   = {...
                                        'delta1', '1, 2',   'mean'; ...
                                        'delta2', '2, 4',   'mean'; ...
                                        'theta',  '4, 8',   'mean'; ...
                                        'alpha',  '8, 13',  'mean'; ...
                                        'beta1',  '13, 20', 'mean'; ...
                                        'beta2',  '20, 30', 'mean'; ...
                                        'gamma1', '30, 40', 'mean'};

% === SCOUTS
sProcess.options.scouts.Comment = '';
sProcess.options.scouts.Type    = 'scout';
sProcess.options.scouts.Value   = {};
sProcess.options.scouts.Hidden  = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    % Get type of data
    Comment = [sProcess.Comment, ':'];
    % Get selected scouts
    ScoutsList = sProcess.options.scouts.Value;
    % Get scouts names
    if ~isempty(ScoutsList) && iscell(ScoutsList) && (size(ScoutsList, 2) >= 2) && ~isempty(ScoutsList{1,2}) && iscell(ScoutsList{1,2})
        ScoutsNames = ScoutsList{1,2};
    elseif ~isempty(ScoutsList) && isstruct(ScoutsList)
        ScoutsNames = {ScoutsList.Label};
    else
        ScoutsNames = [];
    end
    % Format comment
    if isempty(ScoutsNames)
        Comment = [Comment, ' [no selection]'];
    else
        if (length(ScoutsNames) > 15)
            Comment = [Comment, ' [', num2str(length(ScoutsNames)), ' scouts]'];
        else
            for i = 1:length(ScoutsNames)
                Comment = [Comment, ' ', ScoutsNames{i}];
            end
        end
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

id_method = sProcess.options.method.Value{1}; % acquisisco il metodo
method    = sProcess.options.method.Value{2}{id_method};

id_coupling = sProcess.options.coupling.Value{1}; % acquisisco il tipo di coupling 
couplings   = {'<none>', 'tot', 'intraHemi', 'interHemi', 'intraPreRol', 'intraPostRol', 'PrevsPostRol'};
coupling    = couplings{id_coupling};
clear couplings;

foldername = strtrim(sProcess.options.folder.Value);
subfoldername = strtrim(sProcess.options.subfolder.Value);

ntrials     = sProcess.options.boots.Value{1};
perc        = sProcess.options.perctl.Value{1};
nameCat     = strsplit(sProcess.options.nameCat.Value,',');
dataNames   = {'_obs','_sur'};

rng('shuffle');

pBar  = bst_progress('start','Evaluating Source Coupling','Evaluating Source Coupling',0,numel(sInputs));
scrsz = get(0,'ScreenSize'); scrsz(2) = scrsz(3)/2; scrsz = scrsz.*[1 1 1/2 1/2];

subjectNames = {sInputs.SubjectName};
subjectNames = unique(subjectNames);
nSubjects = numel(subjectNames);

tmpcell = cell([nSubjects,length(dataNames)]);
tmpcellDx = cell([nSubjects,length(dataNames)]);
tmpcellSx = cell([nSubjects,length(dataNames)]);

for subjectIdx = 1:nSubjects
    
    for dataTypeIdx = 1:numel(dataNames)
        
        subjectMask = ~cellfun(@isempty,regexp({sInputs.FileName},unique(subjectNames(subjectIdx)),'match'));
        typeMask= ~cellfun(@isempty,regexp({sInputs.Comment},dataNames(dataTypeIdx),'match'));
        fileIdx = find(subjectMask & typeMask);
        
        % === LOAD SURFACE ===
        sSubject = bst_get('Subject', sInputs(subjectIdx).SubjectFile);
        % Get default cortex surface 
        SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
        % Load surface
        sSurf = in_tess_bst(SurfaceFile);
        if isempty(sSurf) || ~isfield(sSurf, 'Atlas')
            bst_report('Error', sProcess, sInputs(iInput), ['Invalid surface file: ' SurfaceFile]);
        end
        
        % Get scouts
        AtlasList = sProcess.options.scouts.Value;
        % Convert from older structure (keep for backward compatibility)
        if isstruct(AtlasList) && ~isempty(AtlasList)
            AtlasList = {'User scouts', {AtlasList.Label}};
        end
        % No scouts selected: exit
        if isempty(AtlasList) || ~iscell(AtlasList) || (size(AtlasList,2) < 2) || isempty(AtlasList{1,2})
            bst_report('Error', sProcess, [], 'No scout selected.');
            return;
        end

        % Get the name of the first atlas obtained from the GUI (ignore the
        % others since there should be exactly one selected) and extract
        % its scouts.
        AtlasName = AtlasList{1,1};
        iAtlasSurf = find(strcmpi(AtlasName, {sSurf.Atlas.Name}));
        CompleteAtlas = {sSurf.Atlas(iAtlasSurf).Scouts.Label}; %#ok<FNDSB>
        
        % Loop on the scouts selected for this atlas
        idx_good   = false(length(CompleteAtlas),1);
        for iScout = 1:length(AtlasList{1,2})
            % Get scout name and locate it in the complete atlas
            ScoutName = AtlasList{1,2}{iScout};
            idx_good( find(strcmpi(ScoutName,CompleteAtlas), 1) ) = true;
        end
        
        %sMat = in_bst_results(sInputs(iInput).FileName, 0);
        % read one file at time
        data = in_bst_timefreq(sInputs(1,fileIdx).FileName); %#ok<FNDSB>
        nFreqs = size(data.Freqs,1); % numero di frequenze
        
        if any(strcmpi(method, {'corr','cohere','plv','plvt','aec','wpli'})) % metric is symmetric
            data.TF = process_compress_sym('Expand', data.TF, (sqrt(8*size(data.TF,1)+1)-1)/2);
        end
        
        nSources = sqrt(size(data.TF,1)); % number subject's channels
        
        data_RowNames = data.RowNames(idx_good);
        tmp = reshape(data.TF, nSources, nSources, nFreqs);
      
        switch coupling
            case 'tot'
                tmp = tmp(idx_good,idx_good,:);
                tmpcell{subjectIdx,dataTypeIdx} = tmp;
                
            case 'intraHemi'
                id_mask = compute_mask(coupling, sScouts); 
                idSx = ismember(id_mask(:,1),find(~idx_good)); % guardo se ci sono canali brutti 
                idDx = ismember(id_mask(:,2),find(~idx_good));
                id_mask(idSx | idDx,:) = []; % cerco gli ID dei bad channel per eliminarli in entrambe le colonne delle matrici

                tmpSx = tmp(id_mask(:,1),id_mask(:,1),:);            
                tmpDx = tmp(id_mask(:,2),id_mask(:,2),:);
                tmpcellSx{subjectIdx,dataTypeIdx} = tmpSx; % creo la matrice per emisfero sinistro
                tmpcellDx{subjectIdx,dataTypeIdx} = tmpDx; % creo la matrice per emisfero destro
                
            case 'interHemi'
                id_mask = compute_mask(coupling, sScouts); 
                idSx   = ismember(id_mask(:,1),find(~idx_good)); % guardo se ci sono canali brutti 
                idDx   = ismember(id_mask(:,2),find(~idx_good));
                id_mask(idSx | idDx,:) = [];
                tmp = tmp(id_mask(:,1), id_mask(:,2),:);
                tmpcell{subjectIdx,dataTypeIdx} = tmp;
                
            case 'intraPreRol'
                idmask = compute_mask(coupling); 
                id_mask = setdiff(idmask, find(~idx_good));
                tmp = tmp(id_mask, id_mask,:);
                tmpcell{subjectIdx,dataTypeIdx} = tmp;
                
            case 'intraPostRol'
                idmask = compute_mask(coupling);
                id_mask = setdiff(idmask, find(~idx_good));
                tmp = tmp(id_mask, id_mask,:);
                tmpcell{subjectIdx,dataTypeIdx} = tmp;
                
            case 'PrevsPostRol'
                id_mask = compute_mask(coupling);
                idPre  = ismember(id_mask(:,1),find(~idx_good)); % guardo se ci sono canali brutti 
                idPost = ismember(id_mask(:,2),find(~idx_good));
                id_mask(idPre | idPost,:) = [];
                tmp = tmp(id_mask(:,1), id_mask(:,2),:);
                tmpcell{subjectIdx,dataTypeIdx} = tmp;
        end
    end
end

if contains('intraHemi', coupling)
    obsDx = tmpcellDx(~cellfun('isempty', tmpcellDx(:,1)),1);
    surDx = tmpcellDx(~cellfun('isempty', tmpcellDx(:,2)),2);
   
    [~, ~, meanCatObsDx, meanCatSurDx, percCatObsDx, percCatSurDx] = compute_Coupling(obsDx, surDx, method, coupling, nameCat, ntrials, perc); 
    
    obsSx = tmpcellSx(~cellfun('isempty', tmpcellSx(:,1)),1);
    surSx = tmpcellSx(~cellfun('isempty', tmpcellSx(:,2)),2);
    
    [~, ~, meanCatObsSx, meanCatSurSx, percCatObsSx, percCatSurSx] = compute_Coupling(obsSx, surSx, method, coupling, nameCat, ntrials, perc);
    
    ncat = length(nameCat);
    meanCatObs  = zeros(ncat,nFreqs);
    meanCatSur  = zeros(ncat,nFreqs);
    percCatObs  = zeros(length(perc), nFreqs, ncat);
    percCatSur  = zeros(length(perc), nFreqs, ncat);
    darkblue  = [0.00, 0.00, 0.55];
    darkgreen = [0.00, 0.39, 0.00];
    darkred   = [0.55, 0.00, 0.00];
    colorsObs   = {darkred, darkgreen, darkblue};
    
    figure('Position',scrsz);
    ax = gca;

    for c = 1:ncat
        
        meanCatObs(c,:)   = (meanCatObsDx(c,:) + meanCatObsSx(c,:))./2;
        meanCatSur(c,:)   = (meanCatSurDx(c,:) + meanCatSurSx(c,:))./2;
        
        for pc = 1:length(perc)
            percCatObs(pc,:,c) = (percCatObsDx(pc,:,c) + percCatObsSx(pc,:,c))./2;
            percCatSur(pc,:,c) = (percCatSurDx(pc,:,c) + percCatSurSx(pc,:,c))./2;
        end
        
        plotfillbetween(ax, percCatObs(:,:,c), meanCatObs(c,:), colorsObs{c}, '-'); hold on;
        ph = plot(ax, meanCatSur(c,:));
        ph.LineStyle = '--';
        ph.Color     = colorsObs{c};
       % plotfillbetween(ax, percCatSur(:,:,c), meanCatSur(c,:), colorsObs{c}, '--'); hold on;      
          
    end
    
    title(ax, [ upper(method) ' (' coupling ') final']);
    %xlim([1 40])
    cts = length(nameCat);
    lgnd = cell(2*cts,1);
    for ct = 1:length(nameCat)
        lgnd{(ct-1)*2+1} = [nameCat{ct}  ' obs'];
        lgnd{      ct*2} = [nameCat{ct} ' sur'];
    end
    
    legend(ax,lgnd);
    hold off;
    OutputFiles = [];
else
    
    obs = tmpcell(~cellfun('isempty', tmpcell(:,1)),1);
    sur = tmpcell(~cellfun('isempty', tmpcell(:,2)),2);
    
    [listobs, listsur] = compute_Coupling(obs, sur, method, coupling, nameCat, ntrials, perc);
    
    if strcmp(coupling, 'tot')
        % save lisobs e listsur in OutputFiles
        if isempty(bst_get('Subject', foldername, 1))
            db_add_subject(foldername);
        end
        
        iStudy                = db_add_condition(foldername,subfoldername);
        
        DataMatObs.TF		  = listobs;
        DataMatObs.Comment    = strcat(method,' listobs');
        DataMatObs.Time		  = [min(data.Time) max(data.Time)];
        DataMatObs.Freqs      = data.Freqs;
        DataMatObs.DataType   = 'data';
        DataMatObs.Method     = method;
        DataMatObs.DataFile   = [];
        DataMatObs.Measure    = 'other';
        DataMatObs.RefRowNames= data_RowNames;
        DataMatObs.RowNames   = data_RowNames;
        DataMatObs.nAvg       = length(sInputs);
        
        DataMatSur.TF		  = listsur;
        DataMatSur.Comment    = strcat(method,' listsur');
        DataMatSur.Time		  = [min(data.Time) max(data.Time)];
        DataMatSur.Freqs      = data.Freqs;
        DataMatSur.DataType   = 'data';
        DataMatSur.Method     = method;
        DataMatSur.DataFile   = [];
        DataMatSur.Measure    = 'other';
        DataMatSur.RefRowNames= data_RowNames;
        DataMatSur.RowNames   = data_RowNames;
        DataMatSur.nAvg       = length(sInputs);
        
        OutputFiles{1} = bst_process('GetNewFilename', ...
            fileparts([foldername, '/', subfoldername, '/']), ['timefreq_connectn_', method, '_k_obs']);
        OutputFiles{2} = bst_process('GetNewFilename', ...
            fileparts([foldername, '/', subfoldername, '/']), ['timefreq_connectn_', method, '_k_sur']);
        
        OutputFiles = cleanConflicts(OutputFiles);
        save(OutputFiles{1}, '-struct', 'DataMatObs');
        save(OutputFiles{2}, '-struct', 'DataMatSur');
        
        db_add_data(iStudy, OutputFiles{1}, DataMatObs);
        db_add_data(iStudy, OutputFiles{2}, DataMatSur);
    else    
        OutputFiles = [];
    end
end
end

%% SUPPORT FUNCTIONS

function OutputFiles = cleanConflicts(OutputFiles)

    if OutputFiles{1} == OutputFiles{2}
        strTmp = OutputFiles{2};
        strTmp(end-8:end) = strcat('_',regexprep(strcat(num2str(randi(9,4,1)')),'\s+',''),'.mat');
        OutputFiles{2} = strTmp;
    end   
end


function mask = compute_mask(coupling, sScouts)

nScouts = length(sScouts);
idDx = [];
idSx = [];

for k = 1:nScouts
    if sScouts(k).Region(1) == 'R'
        idDx = [idDx, k];
    elseif sScouts(k).Region(1) == 'L'
        idSx = [idSx, k];
    else
        error('Atlas contains scouts in ''U'' hemisphere (not ''R'' nor ''L'').');
    end
end

%chansDx   = {'Fp2', 'F4',  'C4',  'P4', 'O2', 'F8', 'T4', 'T6'};
%idDx      = [    2,    4,     6,     8,   10,   12,   14,   16];

%chansSx   = {'Fp1', 'F3',  'C3',  'P3', 'O1', 'F7', 'T3', 'T5'};
%idSx      = [    1,    3,     5,     7,    9,   11,   13,   15];

%chansCx   = {'Fz',  'Cz',  'Pz'};
%idCx      = [  17,    18,    19];

%chansPre  = {'Fp1', 'Fp2', 'F3',  'F4', 'F7', 'F8'};
idPre     = [    1,     2,    3,     4,   11,   12];

%chansPost = {'P3',  'P4',  'O1',  'O2', 'T5', 'T6'};
idPost    = [   7,     8,     9,    10,   15,   16];

%chansRol  = {'C3',  'C4',  'T3',  'T4'};
%idRol     = [   5,     6,    13,    14];

switch coupling
    case {'intraHemi', 'interHemi'}
          mask(:,1) = idSx; 
          mask(:,2) = idDx;
    
    case 'intraPreRol'
          mask(:,1) = idPre;
    
    case 'intraPostRol'
         mask(:,1) = idPost;
    
    case 'PrevsPostRol'
        mask(:,1) = idPre;
        mask(:,2) = idPost;
end
end

function [listobs, listsur, meanCatObs,  meanCatSur, percCatObs, percCatSur] = compute_Coupling(obs, sur, method, coupling, nameCat, ntrials, perc)

[listobs,listsur] = compute_metric(obs, sur, method);

%% SETUP FOR BOOTSTRAPPING

ncat   = length(nameCat);
nSubjs = size(listobs,1);
assert(mod(nSubjs,ncat) == 0);
np     = nSubjs/ncat;
nFreqs = size(listobs{1},3);

meanCatObs  = zeros(ncat,nFreqs);
meanCatSur  = zeros(ncat,nFreqs);
percCatObs  = zeros(length(perc), nFreqs, ncat);
percCatSur  = zeros(length(perc), nFreqs, ncat);

darkblue  = [0.00, 0.00, 0.55];
darkgreen = [0.00, 0.39, 0.00];
darkred   = [0.55, 0.00, 0.00];
colorsObs = {darkred, darkgreen, darkblue};
scrsz     = get(0,'ScreenSize'); scrsz(2) = scrsz(3)/2; scrsz = scrsz.*[1 1 1/2 1/2];

figure('Position',scrsz);
ax = gca;

%% BOOTSTRAPPING

for c = 1:ncat
    
    subj = (c-1)*np+1:c*np;
    
    bootMeanObs = bootstrap_subjects(listobs(subj), ntrials);
    bootMeanSur = bootstrap_subjects(listsur(subj), ntrials);
    
    meanCatObs(c,:)   = compute_mean(listobs(subj));
    meanCatSur(c,:)   = compute_mean(listsur(subj));
    
    percCatObs(:,:,c) = prctile(bootMeanObs, perc);
    percCatSur(:,:,c) = prctile(bootMeanSur, perc);
    
    plotfillbetween(ax, percCatObs(:,:,c), meanCatObs(c,:), colorsObs{c}, '-'); hold on;
    ph = plot(ax, meanCatSur(c,:));
    ph.LineStyle = '--';
    ph.Color     = colorsObs{c};
end

%xlim([1 40]);
if strcmp(method, 'plv')
    title(ax, ['i' upper(method) ' (' coupling ')']);
else
    title(ax, [upper(method) ' (' coupling ')']);
end
xlabel('Frequency bands');

cts = length(nameCat);
lgnd = cell(2*cts,1);
for ct = 1:length(nameCat)
    lgnd{(ct-1)*2+1} = [nameCat{ct}  ' obs'];
    lgnd{      ct*2} = [nameCat{ct} ' sur']; 
end

legend(ax,lgnd);
hold off;

end

function meanCat = compute_mean(lst)

nSubj    = length(lst);
nFreq    = size(lst{1},3);
meansubj = zeros(nSubj,nFreq);
    
for k = 1:nSubj
    meansubj(k,:) = mean( reshape(lst{k}, [], nFreq), 1 );
end

meanCat = mean(meansubj);

end

function [ListObs,ListSur] = compute_metric(o,s,method)

    nSubjs  = length(o);
    ListObs = cell(nSubjs,1);
    ListSur = cell(nSubjs,1);
    
    for k = 1:nSubjs

        switch method
            case 'plv'
                f = @(x) abs(imag(x));
            case 'wpli'
                f = @abs;
            case 'oCC'
                f = @abs;

        end

        obs = o{k};
        sur = s{k};
        
        for ff = 1:size(obs,3)
            obs(:,:,ff) = obs(:,:,ff) - diag(diag(obs(:,:,ff)));
            sur(:,:,ff) = sur(:,:,ff) - diag(diag(sur(:,:,ff)));
        end
        
        ListObs{k,1} = f(obs);
        ListSur{k,1} = f(sur);

    end
end


function bootcomp = bootstrap_subjects(listsubj,ntrials)

    nSubjs = size(listsubj,1);
    nFreqs = size(listsubj{1},3);
    meanchan = zeros(nSubjs,nFreqs);
    bootcomp = zeros(ntrials,nFreqs);

    for i = 1:ntrials
        idxsubj = randi([1 nSubjs],1,nSubjs);
        for j = 1:nSubjs
            nChans = size(listsubj{idxsubj(j)},1);
            subj = reshape(listsubj{idxsubj(j)}, nChans*nChans, nFreqs);
            meanchan(j,:) = mean(subj,1);
        end
        bootcomp(i,:) = mean(meanchan);
    end

end

function [ph] = plotfillbetween(ax, percCat, meanCat, color, linest)
% PLOTFILLBETWEEN Plots bootstrapped mean for two lists of patients
%
%  Obs is a nChans^2 x nFreqs x nSubjs
%
axes(ax);
nFreqs = size(meanCat,2);

x = 1:nFreqs;
y1 = percCat(1,:,1);
y2 = percCat(2,:,1);

ph = plot(x, meanCat(1,:), 'LineStyle', linest);
ph.Color = color;

hold on;

fh = fill([x fliplr(x)], [y1 fliplr(y2)], color, 'HandleVisibility', 'off');

fh.EdgeColor = 'none';
fh.FaceAlpha = 0.15;

hold off;
end



%% ===== GET SCOUTS =====
% USAGE:  [sScouts, sSurf, iSurf] = GetScouts(SurfaceFile)
%         [sScouts, sSurf, iSurf] = GetScouts(iScouts)
%         [sScouts, sSurf, iSurf] = GetScouts()
function [sScouts, sSurf, iSurf] = GetTheScouts(SurfaceFile)
    global GlobalData;
    sScouts = [];
    sSurf   = [];
    iSurf   = [];
    % Parse input
    if (nargin < 1) || isempty(SurfaceFile)
        SurfaceFile = GlobalData.CurrentScoutsSurface;
        iScouts = [];
    elseif ischar(SurfaceFile)
        iScouts = [];
    else
        iScouts = SurfaceFile;
        SurfaceFile = GlobalData.CurrentScoutsSurface;
    end
    % If no surface is defined : do nothing
    if isempty(SurfaceFile)
        return
    end
    % Get loaded surface
    [sSurf, iSurf] = bst_memory('LoadSurface', SurfaceFile);
    % Get the selected scouts
    if (length(sSurf) == 1) && ~isempty(sSurf.Atlas) && ~isempty(sSurf.iAtlas)
        sScouts = sSurf.Atlas(sSurf.iAtlas).Scouts;
        % Select only the required scouts
        if (any(iScouts) > length(sScouts))
            error('Invalid scout indice.');
    	elseif ~isempty(iScouts)
            if any(iScouts > length(sScouts))
                disp('Error: Invalid indices');
            else
                sScouts = sScouts(iScouts);
            end
        end
    end
end


