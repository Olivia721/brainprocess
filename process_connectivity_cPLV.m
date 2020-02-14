function varargout = process_connectivity_cPLV( varargin )

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
sProcess.Comment     = 'Phase Synchrony a.f.o. Frequency';
sProcess.FileTag     = '__';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Connectivity';
sProcess.Index       = 901;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'timefreq'};
sProcess.OutputTypes = {'timefreq'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(~, sInputs) %#ok<DEFNU>

profile on;
% init output PLV
DataMatObserved = db_template('timefreqmat');
DataMatSurrogate = db_template('timefreqmat');
method = 'plv';

rng('shuffle');

bst_progress('start','Evaluating Phase Sync','Evaluating Phase Sync',0,numel(sInputs));

for fileIdx = 1:numel(sInputs)
    % read one file at time
    data = in_bst_timefreq(sInputs(fileIdx).FileName);
    channels = in_bst_channel(sInputs(fileIdx).ChannelFile);
    % TODO fix with option
    iChannels = channel_find(channels.Channel, 'EEG');
    
    [nChans, ~, nFreqs] = size(data.TF);
    % init output
    nEdges   = nChans*nChans;
    edgeData = zeros(nEdges,2,nFreqs);
    
    if max(iChannels) == nChans
        % for each frequency ...
        for freqIdx = 1:nFreqs
            edgeData(:,:,freqIdx) = computePhaseMetric(data.TF(iChannels,:,freqIdx), method);
        end
    else
        warning('SKIPPING FOR LOOP!!! \n');
    end
    
    % ===== SAVE THE RESULTS =====
    % % Get the output study (pick the one from the first file)
    DataMatObserved.TF			= squeeze(edgeData(:,1,:));
    
    parentStruct                = bst_process('GetInputStruct',data.DataFile);
    
    
    %iStudy                      = sInputs(fileIdx).iStudy;
    DataMatObserved.Comment     = strcat(method,' observed');
    %DataMatObserved.ChannelFlag = data.ChannelFlag;% List of good/bad channels (1=good, -1=bad)
    
    DataMatObserved.Time		= [min(data.Time) max(data.Time)];
    DataMatObserved.Freqs       = data.Freqs;
    DataMatObserved.DataType    = 'data';
    DataMatObserved.Method      = method;
    DataMatObserved.DataFile    = parentStruct.FileName;
    DataMatObserved.Measure     = 'other';
    DataMatObserved.RefRowNames = [{channels.Channel.Name}];
    DataMatObserved.RowNames    = [{channels.Channel.Name}];
    DataMatObserved.nAvg        = length(sInputs); % Number of epochs that were averaged to get this file
    
    DataMatSurrogate.TF			= squeeze(edgeData(:,2,:));
    
    iStudy                      = sInputs(fileIdx).iStudy;
    DataMatSurrogate.Comment    = strcat(method,' surrogate');
    %DataMatSurrogate.ChannelFlag = data.ChannelFlag;% List of good/bad channels (1=good, -1=bad)
    
    DataMatSurrogate.Time		= [min(data.Time) max(data.Time)];
    DataMatSurrogate.Freqs      = data.Freqs;
    DataMatSurrogate.DataType   = 'data';
    DataMatSurrogate.Method     = method;
    DataMatSurrogate.DataFile   = parentStruct.FileName;
    DataMatSurrogate.Measure    = 'other';
    DataMatSurrogate.RefRowNames= [{channels.Channel.Name}];
    DataMatSurrogate.RowNames   = [{channels.Channel.Name}];
    DataMatSurrogate.nAvg       = length(sInputs); % Number of epochs that were averaged to get this file
    
    % Create a default output filename
    OutputFiles{1} = bst_process('GetNewFilename', ...
        fileparts(sInputs(fileIdx).FileName), ['timefreq_connectn_', method, '_obs']);
    OutputFiles{2} = bst_process('GetNewFilename', ...
        fileparts(sInputs(fileIdx).FileName), ['timefreq_connectn_', method, '_sur']);
    
    % required one saving two files in such short time 
    % to avoid conflicts in filenames (BS uses clock to gen. unique
    % filenames)
    OutputFiles = cleanConflicts(OutputFiles);
    
    % Save on disk
    save(OutputFiles{1}, '-struct', 'DataMatObserved');
    save(OutputFiles{2}, '-struct', 'DataMatSurrogate');
    
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, DataMatObserved);
    db_add_data(iStudy, OutputFiles{2}, DataMatSurrogate);
    
    % update progress bar accordingly
    bst_progress('inc',1);
end

profile viewer;
end

%% SUPPORT FUNCTIONS
function OutputFiles = cleanConflicts(OutputFiles)

    if OutputFiles{1} == OutputFiles{2}
        strTmp = OutputFiles{2};
        strTmp(end-8:end) = strcat('_',regexprep(strcat(num2str(randi(9,4,1)')),'\s+',''),'.mat');
        OutputFiles{2} = strTmp;
    end   
end

function [edgeData] = computePhaseMetric(data,metric)
% Description
% 			pTE can but requires an ad hoc function to estimate MI
%	[PTE,PLV, F] = computePhaseMetric(data,channelNames)

switch metric
    case 'plv'
        edgeData = computePlv_fast(data);
    case 'wpli'
        edgeData = computewPLI(data);
    case 'pTE'
        edgeData = [];
    case 'oCC'
        edgeData = computeoCC(data);
    otherwise
        bst_error(sprintf('metric %s unsupported', metric'));
    
    % Phase amplitude coupling??? Amplitude in high freqs. correlates
    % with phase async. for low freqs.
end

end

function [out] = computePlv_fast(signals)

signals = signals ./ abs(signals);

[nChans, nSamples] = size(signals);

plv     = (signals * signals') / nSamples;
plvsup  = triu(plv,1);
plv     = plvsup + plvsup.' + eye(nChans);

I        = randi(nSamples, [nChans,1]);
signals_ = [signals, signals];
signals_ = signals_( bsxfun(@plus, ...
                     nChans*bsxfun(@plus, nSamples-I, 0:(nSamples-1)), ...
                     (1:nChans)'));

plvs   = (signals * signals_') / nSamples;
plvsup = triu(plvs,1);
plvs   = plvsup + plvsup.' + eye(nChans);

out(:,:,1) = plv(:);
out(:,:,2) = plvs(:);

end

% function [out] = computePlv(Xfc)
% % Description
% %	[PLV,PLVSURR] = computePlv(signals,nPermutation)
% 
% % Inst. phase [ch x T]
% 
% Xfc = Xfc./abs(Xfc);
% 
% [nChans,nSamples] = size(Xfc);
% 
% % allocate fullmat plv mat
% cPlv = complex(eye(nChans),zeros(nChans));
% cPlvSurr = complex(eye(nChans),zeros(nChans));
% 
% %compute only lower triangular
% for iCh = 1:nChans-1
%     Xi = Xfc(iCh,:);
%     for jCh = iCh+1:nChans
%         Xj = Xfc(jCh,:);
%         cPlv(iCh,jCh) = (Xi*Xj.')./nSamples;
%         offset = randi(nSamples,1);
%         cPlvSurr(iCh,jCh) = (Xi*circshift(Xj,offset).')./nSamples;
%     end
% end
% 
% 
% %% FASTER CODE? %%
% %parfor linIdx = 2:nEdges+1
% %	[iCh, jCh] = ind2sub([nChans,nChans],linIdx);
% %	cPlv(linIdx) = (Xfc(iCh,:)*Xfc(jCh,:)')./nSamples;
% %  offset = randi(nSamples,1);
% %	cPlvSurr(linIdx) = (Xfc(iCh,:)*circshift(Xfc(jCh,:),offset)')./nSamples;
% %end
% 
% % symmetrize
% cPlv     = cPlv     + triu(cPlv,1).';
% cPlvSurr = cPlvSurr + triu(cPlvSurr,1).';
% 
% % get vectorized mat
% out(:,:,1) = cPlv(:);%(cPlv ~= 0);
% out(:,:,2) = cPlvSurr(:);%(cPlv ~= 0);
% 
% end

function [out] = computewPLI(Xfc)
% Description
%	[PLV,PLVSURR,AMPCORR,AMPCORRSURR] = computewPLI(signals,nPermutation)

% Inst. phase [ch x T]
%Xfc = hilbert(signals')';

[nChans,nSamples] = size(Xfc);

% allocate fullmat plv mat
wPLI     = eye(nChans);
wPLISurr = eye(nChans);

for iCh = 1:nChans-1
    for jCh = iCh+1:nChans
        % real val
        %crossSpectr = (Xfc(iCh,:).*Xfc(jCh,:)'); FABIOs
        crossSpectr = Xfc(iCh,:) .* Xfc(jCh,:);
        X = imag(crossSpectr);
        wPLI(iCh,jCh) = sum(abs(X).* sign(X),2)./ sum(abs(X),2);
        
        % surrogate
        offset = randi(nSamples,1);
        %crossSpectr = (Xfc(iCh,:).*circshift(Xfc(jCh,:),offset)'); FABIO
        crossSpectr = Xfc(iCh,:) .* circshift(Xfc(jCh,:),offset);
        X = imag(crossSpectr);
        wPLISurr(iCh,jCh) = sum(abs(X).* sign(X),2)./ sum(abs(X),2);
    end
end

% get vectorized mat
out(:,1) = wPLI(:);%(wPLI ~= 0);
out(:,2) = wPLISurr(:);%(wPLI ~= 0);
end

function [out] = computeoCC(signals)
% Description
%	[PLV,PLVSURR,AMPCORR,AMPCORRSURR] = computewPLI(signals,nPermutation)

% Inst. phase [ch x T]
%Xfc = hilbert(signals')';
Xfc = signals;

[nChans,nSamples] = size(Xfc);


% allocate fullmat plv mat
oCC = zeros(nChans);
oCCSurr = zeros(nChans);

for iCh = 1:nChans
    
    X    = Xfc(iCh,:);
    absX = abs(X);
    avgX = mean(absX);
    stdX = std(absX);
    XX   = conj(X) ./ absX;
    
    for jCh = 1:nChans

        % orthogonalize Y to X
        Y    = Xfc(jCh,:);
        absY = abs(Y);
        
        Yabsort = imag(absY .* XX);
        
        avgYao  = mean(Yabsort);
        stdYao  = std(Yabsort);
        
        % compute CC with orthogonalized data
        oCC(iCh,jCh) = (absX * Yabsort.' - nSamples*avgX*avgYao) / (nSamples-1) / stdX / stdYao;
        
        % cut and rotate
        % recompute CC under H0
        offset = randi(nSamples,1);
        oCCSurr(iCh,jCh) = (circshift(absX,offset) * Yabsort.' - nSamples*avgX*avgYao) / (nSamples-1) / stdX / stdYao;
        
    end
end

% get vectorized mat
out(:,:,1) = oCC(:);%(oCC ~= 0);
out(:,:,2) = oCCSurr(:);%(oCCSurr ~= 0);
end


