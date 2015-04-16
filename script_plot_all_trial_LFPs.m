% script to plot individual trial LFPs for a single tetrode-session on a
% single sheet

chDB_directory    = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal data structures';
lfp_root          = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/stop-signal LFPs';
saveDir           = '/Volumes/PublicLeventhal1/dan/stop-signal reanalysis/all_LFP_plots';
[chDB_list, chDB_fnames, ~, ~] = get_chStructs_for_analysis;

pages_per_file = 50;

cp = initChanParams();

locs_to_skip = {'ref','emg','eeglam'};
channels_to_skip = {'e2', 'e3'};
timeWindow = [-1 3];

trialType = 'correctgo';
eventName = 'noseCenterIn';
skipLFPs = 1;
colorList = ['r', 'b'];
linewidth = .5;
highlightwidth = 1.5;
yRange = .3;


for i_chDB = 11 : length(chDB_list)
    
    % first, load the relevant channel DBs, if necessary
    if ~exist(chDB_list{i_chDB}, 'var')
        chDB_file = fullfile(chDB_directory, chDB_fnames{i_chDB});
        disp(['loading ' chDB_file]);
        load( chDB_file );
    end
    
    if i_chDB < 5
        implantID = implantID_from_ratID(chDB_list{i_chDB}(1:3));
        chDB_info = whos( [chDB_list{i_chDB}(1:3) 'Ch*'] );
    else
        implantID = chDB_list{i_chDB}(1:5);
        chDB_info = whos( [implantID 'Ch*'] );
    end
    ch = eval( chDB_info.name );
    
    subjectFolder = fullfile(saveDir, [implantID '_lfpPlots']);
    
    base_saveName = [implantID 'LFPs'];
    base_saveName = fullfile(subjectFolder, base_saveName);
    
    lfp_directory = fullfile(lfp_root, [implantID '_LFPs']);
    cd(lfp_directory);

    % figure out where the breaks in the plot files will be
    numValidPlots = 0;
    for iCh = 1 : length(ch)

        c = ch{iCh};

        if any(strcmpi(c.location.subclass, locs_to_skip))
            continue;
        end
        if any(strcmpi(c.tetrode.name, channels_to_skip))
            continue;
        end
        
        sessionDate = c.date;
        if length(sessionDate) > 8
            sessionDate = datestr(sessionDate, 'yyyymmdd');
        end
        lfp_fileinfo = dir(['*_' sessionDate '*.hsdf']);
        if isempty(lfp_fileinfo); continue; end

        lfpFn = lfp_fileinfo.name;
%         if lfpFn(1) == '\'
%             lfpFn = PCfn2macfn( lfpFn );
%         end
%         [~, lfpFn, ext, ~] = fileparts(lfpFn);
%         lfpFn = [lfpFn ext 'f'];
% 
%         lfpFn = which(lfpFn);
%         if isempty(lfpFn)
%             disp(['No lfp file for channel ' c.name]);
%             continue;
%         end

        repWire = getRepWire(c);
        if ~repWire
            disp(['No representative wire for channel ' c.name]);
            continue;
        end
        
        numValidPlots = numValidPlots + 1;
        
        if rem(numValidPlots, pages_per_file) == 1
            fileNum = ceil(numValidPlots / pages_per_file);
            saveName = sprintf([base_saveName '_%03d.pdf'], fileNum);
            
            if ~exist(saveName, 'file')   % this file hasn't been saved yet
                start_ch_idx = iCh;
                break;
            end
        end
    end
    
    numValidPlots = numValidPlots - 1;
    for iCh = start_ch_idx : length(ch)

        c = ch{iCh};

        if any(strcmpi(c.location.subclass, locs_to_skip))
            continue;
        end
        if any(strcmpi(c.tetrode.name, channels_to_skip))
            continue;
        end
        
        sessionDate = c.date;
        if length(sessionDate) > 8
            sessionDate = datestr(sessionDate, 'yyyymmdd');
        end
        lfp_fileinfo = dir(['*_' sessionDate '*.hsdf']);
        if isempty(lfp_fileinfo); continue; end

        lfpFn = lfp_fileinfo.name;
%         if lfpFn(1) == '\'
%             lfpFn = PCfn2macfn( lfpFn );
%         end
%         [~, lfpFn, ext, ~] = fileparts(lfpFn);
%         lfpFn = [lfpFn ext 'f'];
% 
%         lfpFn = which(lfpFn);
%         if isempty(lfpFn)
%             disp(['No lfp file for channel ' c.name]);
%             continue;
%         end

        repWire = getRepWire(c);
        if ~repWire
            disp(['No representative wire for channel ' c.name]);
            continue;
        end

    %     h_fig = figure;
        lfpHeader = getHSDHeader(lfpFn);
        Fs = lfpFs(lfpHeader);

        trialEventParams = getTrialEventParams(trialType);
        events = DKLgetTrialEvents2( c, trialEventParams );
        ts = events.(eventName);

        lfpSnips = extractERPsfromFile( ts, lfpFn, repWire, 'timeWindow', ...
            timeWindow );

        lfpSnips2 = lfpSnips([1 : skipLFPs : size(lfpSnips, 1)], :);

        trialNum = 0;
        for iTrial = 1 : length(c.trials)
    %         iTrial
            if isnan(c.trials(iTrial).correct)
                continue;
            end
            if ~c.trials(iTrial).correct
                continue;
            end
            trialNum = trialNum + 1;

            if c.trials(iTrial).trialType == 2    % go/nogo
                trialEvents{trialNum}(1) = c.trials(iTrial).timestamps.whiteNoise(1) - ...
                    c.trials(iTrial).timestamps.noseCenterIn(1);
            else
                trialEvents{trialNum}(1) = c.trials(iTrial).timestamps.tone(1) - ...
                    c.trials(iTrial).timestamps.noseCenterIn(1);
            end

            switch c.task 
                case 1,
                    trialEvents{trialNum}(2) = c.trials(iTrial).timestamps.whiteNoise(1) - ...
                        c.trials(iTrial).timestamps.noseCenterIn(1);
                    trialEvents{trialNum}(3) = c.trials(iTrial).timestamps.noseCenterOut(1) - ...
                        c.trials(iTrial).timestamps.noseCenterIn(1);
                    trialEvents{trialNum}(4) = c.trials(iTrial).timestamps.noseSideIn(1) - ...
                        c.trials(iTrial).timestamps.noseCenterIn(1);
                    if ~isempty(c.trials(iTrial).timestamps.foodHopperClick)
                        trialEvents{trialNum}(5) = c.trials(iTrial).timestamps.foodHopperClick(1) - ...
                            c.trials(iTrial).timestamps.noseCenterIn(1);
                    end
                case 2,
                    trialEvents{trialNum}(2) = c.trials(iTrial).timestamps.noseCenterOut(1) - ...
                        c.trials(iTrial).timestamps.noseCenterIn(1);
                    trialEvents{trialNum}(3) = c.trials(iTrial).timestamps.noseSideIn(1) - ...
                        c.trials(iTrial).timestamps.noseCenterIn(1);
                    if ~isempty(c.trials(iTrial).timestamps.foodHopperClick)
                        trialEvents{trialNum}(4) = c.trials(iTrial).timestamps.foodHopperClick(1) - ...
                            c.trials(iTrial).timestamps.noseCenterIn(1);
                    end
                case 3,
                    trialEvents{trialNum}(2) = c.trials(iTrial).timestamps.noseCenterOut(1) - ...
                        c.trials(iTrial).timestamps.noseCenterIn(1);
                    if isfield(c.trials(iTrial).timestamps, 'noseSideIn')
                        trialEvents{trialNum}(3) = c.trials(iTrial).timestamps.noseSideIn(1) - ...
                            c.trials(iTrial).timestamps.noseCenterIn(1);
                    end
                    if ~isempty(c.trials(iTrial).timestamps.foodHopperClick)
                        trialEvents{trialNum}(4) = c.trials(iTrial).timestamps.foodHopperClick(1) - ...
                            c.trials(iTrial).timestamps.noseCenterIn(1);
                    end
                case 4,
                    trialEvents{trialNum}(2) = c.trials(iTrial).timestamps.noseCenterOut(1) - ...
                        c.trials(iTrial).timestamps.noseCenterIn(1);
                    if isfield(c.trials(iTrial).timestamps, 'noseSideIn')
                        trialEvents{trialNum}(3) = c.trials(iTrial).timestamps.noseSideIn(1) - ...
                            c.trials(iTrial).timestamps.noseCenterIn(1);
                    end
                    if ~isempty(c.trials(iTrial).timestamps.foodHopperClick)
                        trialEvents{trialNum}(4) = c.trials(iTrial).timestamps.foodHopperClick(1) - ...
                            c.trials(iTrial).timestamps.noseCenterIn(1);
                    end
            end
        end


        numTrials = size(lfpSnips2, 1);

        h_fig = figure('unit','inches','position',[1 1 8.5 11]);
        h_header = subplot('position', [.01 .9 .98 .09]);
        h_plotAxes = subplot('position', [.01 .01 .98 .88]);

        [h_plotAxes] = plotEventTriggeredSamples2(lfpSnips2, timeWindow, 'color', 'k', ...[0.3 0.3 0.3], ...
            'colorlist',colorList,'linewidth',linewidth,'eventcolor',[.5 .5 .5], ...
            'yrange', yRange * skipLFPs, 'axes', h_plotAxes, 'trialEvents', trialEvents);

        set(gca,'xtick',timeWindow,'xticklabel',{num2str(timeWindow(1)) num2str(timeWindow(2))},'yticklabel','','ylim',[(-numTrials*skipLFPs*yRange) 1.75],'xtick',[],'ytick',[], ...
            'xlim',timeWindow, 'fontname', 'arial','fontsize',20);

        set(gcf,'name',[c.name ', ' c.location.name ', ' c.location.subclass]);

        vertLine(0, 'color','k','linewidth',2);

        axes(h_header);

        text(.1,.5,[c.name ', task: ' num2str(c.task) ', correct go, location: ' c.location.name ', subclass: ' c.location.subclass ', channel index: ' num2str(iCh)]);

        numValidPlots = numValidPlots + 1;
        if rem(numValidPlots, pages_per_file) == 1
            fileNum = ceil(numValidPlots / pages_per_file);
            saveName = sprintf([base_saveName '_%03d'], fileNum);
%         end
%         if iCh == numValidPlots
            export_fig(saveName,'-pdf','-q101','-painters')
        else
            export_fig(saveName,'-pdf','-q101','-painters', '-append')
        end

        close(h_fig);

    end

end
