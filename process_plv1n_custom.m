function varargout = process_plv1n_custom( varargin )
% PROCESS_PLV1N: Compute the coherence between all the pairs of signals, in one file.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2019 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
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
% Authors: Francois Tadel, 2012-2014

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Phase locking value (PLV) NxN - custom';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Connectivity';
    sProcess.Index       = 660;
    %sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/Connectivity';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data',     'results',  'matrix'};
    sProcess.OutputTypes = {'timefreq', 'timefreq', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;

    % === CONNECT INPUT
    sProcess = process_corr1n('DefineConnectOptions', sProcess, 1);
    % === FREQ BANDS
    sProcess.options.label2.Comment = '<BR><U><B>Estimator options</B></U>:';
    sProcess.options.label2.Type    = 'label';
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
    % === Mirror
    sProcess.options.mirror.Comment = 'Mirror signal before filtering (not recommended)';
    sProcess.options.mirror.Type    = 'checkbox';
    sProcess.options.mirror.Value   = 0;
    % === KEEP TIME
    sProcess.options.keeptime.Comment = 'Keep time information, and estimate the PLV across trials<BR>(requires the average of many trials)';
    sProcess.options.keeptime.Type    = 'checkbox';
    sProcess.options.keeptime.Value   = 0;
    % === PLV METHOD
    sProcess.options.plvmeasure.Comment = {'None (complex)', 'Magnitude', 'Measure:'};
    sProcess.options.plvmeasure.Type    = 'radio_line';
    sProcess.options.plvmeasure.Value   = 2;
    % === OUTPUT
    sProcess.options.label3.Comment = '<BR><U><B>Output configuration</B></U>:';
    sProcess.options.label3.Type    = 'label';
    % === OUTPUT MODE
    sProcess.options.outputmode.Comment = {'Save individual results (one file per input file)', 'Concatenate input files before processing (one file)', 'Save average connectivity matrix (one file)'};
    sProcess.options.outputmode.Type    = 'radio';
    sProcess.options.outputmode.Value   = 1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputA) %#ok<DEFNU>
    % Input options
    OPTIONS = process_corr1n('GetConnectOptions', sProcess, sInputA);
    if isempty(OPTIONS) 
        OutputFiles = {};
        return
    end

    % Keep time or not: different methods
    if sProcess.options.keeptime.Value
        OPTIONS.Method = 'plvt';
    else
        OPTIONS.Method = 'plv';
    end
    % Filtering bands options
    OPTIONS.Freqs = sProcess.options.freqbands.Value;
    OPTIONS.isMirror = sProcess.options.mirror.Value;
    % PLV measure
    if isfield(sProcess.options, 'plvmeasure') && isfield(sProcess.options.plvmeasure, 'Value') && ~isempty(sProcess.options.plvmeasure.Value) 
        switch (sProcess.options.plvmeasure.Value)
            case 1,  OPTIONS.PlvMeasure = 'none';
            case 2,  OPTIONS.PlvMeasure = 'magnitude';
        end
    else
        OPTIONS.PlvMeasure = 'magnitude';
    end
    % Compute metric
    OutputFiles = my_connectivity({sInputA.FileName}, [], OPTIONS);
end



%% ===== TEST FUNCTION =====
function Test() %#ok<DEFNU>
    % Start a new report
    bst_report('Start');
    % Get test datasets
    sFile = process_simulate_ar('Test');
    % Coherence process
    sTmp = bst_process('CallProcess', 'process_plv1n_custom', sFile, [], ...
        'timewindow',   [], ...    % All the time in input
        'freqbands',    {'delta1', '1, 2', 'mean', 'delta2', '2, 4', 'mean'; 'theta', '4, 8', 'mean'; 'alpha', '8, 13', 'mean'; 'beta1', '13, 20', 'mean'; 'beta2', '20, 30', 'mean'; 'gamma1', '30, 40', 'mean'}, ...
        'mirror',       1, ...
        'keeptime',     0, ...
        'outputmode',   1);        % Save individual results (one file per input file)
    % Snapshot: spectrum
    bst_process('CallProcess', 'process_snapshot', sTmp, [], ...
        'target',       11, ...  % Connectivity matrix (image)
        'modality',     1, 'orient', 1, 'time', 0, 'contact_time', [-40, 110], 'contact_nimage', 16, ...
        'Comment',      [sFile.Comment, ': ' sTmp.Comment]);
    % Save and display report
    ReportFile = bst_report('Save', sFile);
    bst_report('Open', ReportFile);
end




