%% Setup & load data
close all; clear; clc;

requiredFiles = ["SouthAsia_DI.mat", "NS_Result_C8.mat", "cond_prob_results.mat"];
for f = requiredFiles
    if ~isfile(f)
        error('Required file "%s" was not found in the current folder.', f);
    end
end

load("SouthAsia_DI.mat")
load("NS_Result_C8.mat")
load("cond_prob_results.mat")

requiredVars = ["SPI","SMI","SRI","cluster_map","cond_prob_results"];
missingVars = requiredVars(~ismember(requiredVars, string(who)));
if ~isempty(missingVars)
    error('Required variables are missing: %s', strjoin(cellstr(missingVars), ', '));
end

dataMap = struct('SPI', {SPI}, 'SMI', {SMI}, 'SRI', {SRI});

cohDir = fullfile("figures", "wavelet_coherence");
lagDir = fullfile("figures", "wavelet_lag");
pdDir = fullfile("figures", "wavelet_probability_density");
dirList = [cohDir, lagDir, pdDir];
for d = dirList
    if ~exist(d, 'dir')
        mkdir(d);
    end
end

%% Run wavelet analysis per cluster
d_type = {'SPI','SMI'; 'SPI','SRI'; 'SMI','SRI'};
alphabet = string(('a':'z').').';
tmp_yr = 1976:2024;

tic
cause_vector = cell(size(d_type,1), max(cluster_map(:)));
effect_vector = cell(size(d_type,1), max(cluster_map(:)));

for typeN = 1:size(d_type,1)
    fig = figure('Position',[2000 0 1200 1000],'Color','w');
    tiledlayout(4,2,'Padding','compact','TileSpacing','compact');

    for clmN = 1:max(cluster_map(:))
        [latN, lonN] = find(cluster_map == clmN);
        if isempty(cond_prob_results{typeN,1}{clmN,2})
            nexttile;
            text(0.5, 0.5, 'No valid data', 'HorizontalAlignment', 'center');
            axis off
            continue
        end

        scaleN = mode(cond_prob_results{typeN,1}{clmN,2});
        cause_vector{typeN,clmN} = [];
        effect_vector{typeN,clmN} = [];

        for grdN = 1:numel(latN)
            cause_vector{typeN,clmN} = [cause_vector{typeN,clmN} squeeze(dataMap.(d_type{typeN,1}){scaleN,1}(latN(grdN),lonN(grdN),13:end))];
            effect_vector{typeN,clmN} = [effect_vector{typeN,clmN} squeeze(dataMap.(d_type{typeN,2}){1,1}(latN(grdN),lonN(grdN),13:end))];
        end

        cause_vector{typeN,clmN} = mean(cause_vector{typeN,clmN},2,'omitnan');
        effect_vector{typeN,clmN} = mean(effect_vector{typeN,clmN},2,'omitnan');

        nt = nexttile; hold on;
        wcoherence(cause_vector{typeN,clmN}, effect_vector{typeN,clmN}, years(1/12));
        colormap('turbo'); colorbar('off'); hold off;
        xlabel('Time','FontWeight','bold'); ylabel('Period (year)','FontWeight','bold');
        xticks(5:10:45); xlim([1 49]); xticklabels(string(tmp_yr(5:10:45)));
        set(gca,'FontName','Arial','FontSize',14);
        if clmN <= 6
            set(gca,'XLabel',[]);
        end
        if mod(clmN,2) ~= 1
            set(gca,'YLabel',[]);
        end
        title(sprintf('(%s) Cluster %d, %s-%d → %s-1', ...
            alphabet(clmN), clmN, d_type{typeN,1}, scaleN, d_type{typeN,2}), ...
            'FontSize', 18, 'FontWeight', 'bold');

        tmp_h = findobj(nt,'Type','patch');
        set(tmp_h,'FaceAlpha',0.7);
        for i = 1:length(tmp_h)
            v = tmp_h(i).Vertices;
            v_new = mean(v) + 0.7*(v - mean(v));
            tmp_h(i).Vertices = v_new;
        end
    end

    saveas(fig, fullfile(cohDir, sprintf("Type%d_wavelet_coherence.png", typeN)));
end
toc

%% Wavelet-derived lag analyses
tmp_date = datetime([repelem(tmp_yr,1,12); repmat(1:12,1,length(tmp_yr)); ones(1,588)]');

tic
for typeN = 1:size(d_type,1)
    fig1 = figure('Position',[0 0 1200 1000],'Color','w');
    tl1 = tiledlayout(fig1,4,2,'Padding','compact','TileSpacing','compact');

    fig2 = figure('Position',[1000 0 1200 500],'Color','w');
    tl2 = tiledlayout(fig2,2,4,'Padding','compact','TileSpacing','compact');

    fig3 = figure('Position',[2000 0 1200 1000],'Color','w');
    tl3 = tiledlayout(fig3,4,2,'Padding','compact','TileSpacing','compact');

    for clmN = 1:max(cluster_map(:))
        if isempty(cond_prob_results{typeN,1}{clmN,2}) || isempty(cause_vector{typeN,clmN}) || isempty(effect_vector{typeN,clmN})
            nexttile(tl1); text(0.5,0.5,'No valid data','HorizontalAlignment','center'); axis off
            nexttile(tl2); text(0.5,0.5,'No valid data','HorizontalAlignment','center'); axis off
            nexttile(tl3); text(0.5,0.5,'No valid data','HorizontalAlignment','center'); axis off
            continue
        end

        scaleN = mode(cond_prob_results{typeN,1}{clmN,2});
        [wcoh, wcs, f, coi] = wcoherence(cause_vector{typeN,clmN}, effect_vector{typeN,clmN}, years(1/12));
        periods = years(f);
        period_range = [0 5];
        coherence_thr = 0.7;
        period_idx = find(periods >= period_range(1) & periods <= period_range(2));

        wcs_region = wcs(period_idx,:);
        wcoh_region = wcoh(period_idx,:);
        periods_region = periods(period_idx);

        time_vector = 1:length(cause_vector{typeN,clmN});
        n_obs = numel(time_vector);
        midpoint_idx = floor(n_obs/2);
        past_indices = 1:midpoint_idx;
        recent_indices = (midpoint_idx+1):n_obs;

        total_lags = [];
        past_lags = [];
        recent_lags = [];

        for t_idx = 1:n_obs
            if time_vector(t_idx) > years(coi(t_idx))
                for p_idx = 1:length(periods_region)
                    if wcoh_region(p_idx, t_idx) > coherence_thr
                        current_phase = angle(wcs_region(p_idx, t_idx));
                        if current_phase > 0
                            lag = (current_phase / (2*pi)) * periods_region(p_idx);
                            total_lags = [total_lags; lag];
                            if ismember(t_idx, past_indices)
                                past_lags = [past_lags; lag];
                            elseif ismember(t_idx, recent_indices)
                                recent_lags = [recent_lags; lag];
                            end
                        end
                    end
                end
            end
        end

        window_size = 120;
        step_size = 12;
        end_time = length(time_vector) - window_size;
        center_times = [];
        median_lags = [];
        percentile_25 = [];
        percentile_75 = [];

        for t = 1:step_size:end_time
            window_indices = t:(t + window_size - 1);
            window_lags = [];
            for t_idx = window_indices
                if time_vector(t_idx) > years(coi(t_idx))
                    for p_idx = 1:length(periods_region)
                        if wcoh_region(p_idx, t_idx) > coherence_thr
                            current_phase = angle(wcs_region(p_idx, t_idx));
                            if current_phase > 0
                                lag = (current_phase / (2*pi)) * periods_region(p_idx);
                                window_lags = [window_lags; lag];
                            end
                        end
                    end
                end
            end
            center_times = [center_times; t + window_size/2];
            median_lags = [median_lags; median(window_lags, 'omitnan')];
            percentile_25 = [percentile_25; prctile(window_lags, 25)];
            percentile_75 = [percentile_75; prctile(window_lags, 75)];
        end

        median_lag_total = median(total_lags, 'omitnan');

        nexttile(tl1);
        if isempty(total_lags)
            text(0.5, 0.5, 'No valid lag data', 'HorizontalAlignment', 'center');
            axis off
        else
            histogram(total_lags, 'Normalization', 'pdf');
            hold on
            xline(median_lag_total, '--r', sprintf('%.2f year\n(%.2f months)', median_lag_total, median_lag_total*12), 'LineWidth', 2);
            hold off
            legend("","median");
            set(gca,'FontName','Arial','FontSize',14)
            xlabel('Propagation Time / Lag (years)','FontWeight','bold'); xlim([0 0.5])
            ylabel('Probability Density','FontWeight','bold');
            if clmN <= 6, set(gca,'XLabel',[]); end
            if mod(clmN,2) ~= 1, set(gca,'YLabel',[]); end
            title(sprintf('(%s) Cluster %d, %s-%d → %s-1', ...
                alphabet(clmN), clmN, d_type{typeN,1}, scaleN, d_type{typeN,2}), ...
                'FontSize', 18, 'FontWeight', 'bold');
        end

        nexttile(tl2);
        if isempty(past_lags) && isempty(recent_lags)
            text(0.5, 0.5, 'No valid lag data', 'HorizontalAlignment', 'center');
            axis off
        else
            group = [repmat({'Past'}, length(past_lags), 1); repmat({'Recent'}, length(recent_lags), 1)];
            data = [past_lags; recent_lags];
            boxplot(data, group,'Labels',{'Past Period', 'Recent Period'},'Widths',0.15,'Symbol','')
            h = findobj(gca,'Tag','Box');
            for boxN = 1:length(h)
                patch(get(h(boxN),'XData'),get(h(boxN),'YData'),[0.4 0.7 1], 'FaceAlpha',0.5,'EdgeColor','none');
            end
            set(findobj(gca,'Tag','Median'),'Color',[0.2 0.4 1],'LineWidth',2)
            set(findobj(gca,'Tag','Whisker'),'LineWidth',2,'Color',[0.2 0.4 1])
            set(gca,'FontName','Arial','FontSize',12)
            ylim([0 0.4]); ylabel('Time / Lag (years)','FontWeight','bold');
            if mod(clmN,4) ~= 1, set(gca,'YLabel',[]); end
            title(sprintf('(%s) Cluster %d\n%s-%d → %s-1', ...
                alphabet(clmN), clmN, d_type{typeN,1}, scaleN, d_type{typeN,2}), ...
                'FontSize', 18, 'FontWeight', 'bold');
        end

        nexttile(tl3);
        if isempty(center_times)
            text(0.5, 0.5, 'No valid lag data', 'HorizontalAlignment', 'center');
            axis off
        else
            hold on
            fill([tmp_date(center_times); flipud(tmp_date(center_times))], [percentile_25; flipud(percentile_75)], ...
                [0.8 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Interquartile Range');
            plot(tmp_date(center_times), median_lags, 'b-', 'LineWidth', 2, 'DisplayName', 'Median Lag');
            hold off
            ylim([0 0.35]); xlim('tight');
            xticks(tmp_date(center_times(1)):years(12):tmp_date(center_times(end)));
            xtickformat('yyyy-MM');
            set(gca,'FontName','Arial','FontSize',14)
            ylabel('Time / Lag (years)','FontWeight','bold');
            if mod(clmN,2) ~= 1, set(gca,'YLabel',[]); end
            title(sprintf('(%s) Cluster %d, %s-%d → %s-1', ...
                alphabet(clmN), clmN, d_type{typeN,1}, scaleN, d_type{typeN,2}), ...
                'FontSize', 18, 'FontWeight', 'bold');
        end
    end

    saveas(fig1, fullfile(lagDir, sprintf("Type%d_lag_histograms.png", typeN)));
    saveas(fig2, fullfile(lagDir, sprintf("Type%d_lag_boxplots.png", typeN)));
    saveas(fig3, fullfile(lagDir, sprintf("Type%d_lag_timeseries.png", typeN)));
end
toc

%% Probability density per cluster
tic
for typeN = 1:size(d_type,1)
    for clmN = 1:max(cluster_map(:))
        if isempty(cond_prob_results{typeN,1}{clmN,2}) || isempty(cause_vector{typeN,clmN}) || isempty(effect_vector{typeN,clmN})
            continue
        end

        scaleN = mode(cond_prob_results{typeN,1}{clmN,2});
        [wcoh, wcs, f, coi] = wcoherence(cause_vector{typeN,clmN}, effect_vector{typeN,clmN}, years(1/12));
        periods = years(f);
        period_range = [0 5];
        coherence_thr = 0.7;
        period_idx = find(periods >= period_range(1) & periods <= period_range(2));

        wcs_region = wcs(period_idx,:);
        wcoh_region = wcoh(period_idx,:);
        periods_region = periods(period_idx);

        time_vector = 1:length(cause_vector{typeN,clmN});
        n_obs = numel(time_vector);
        midpoint_idx = floor(n_obs/2);
        past_indices = 1:midpoint_idx;
        recent_indices = (midpoint_idx+1):n_obs;

        total_lags = [];
        for t_idx = 1:n_obs
            if time_vector(t_idx) > years(coi(t_idx))
                for p_idx = 1:length(periods_region)
                    if wcoh_region(p_idx, t_idx) > coherence_thr
                        current_phase = angle(wcs_region(p_idx, t_idx));
                        if current_phase > 0
                            lag = (current_phase / (2*pi)) * periods_region(p_idx);
                            total_lags = [total_lags; lag];
                        end
                    end
                end
            end
        end

        if isempty(total_lags)
            continue
        end

        median_lag_total = median(total_lags, 'omitnan');
        fig = figure('Position',[2000 500 600 250],'Color','w');
        histogram(total_lags, 'Normalization', 'pdf');
        hold on
        xline(median_lag_total, '--r', sprintf('%.2f year\n(%.2f months)', median_lag_total, median_lag_total*12), 'LineWidth', 2);
        hold off
        legend("","median");
        set(gca,'FontName','Arial','FontSize',14)
        xlabel('Propagation Time / Lag (years)','FontWeight','bold'); xlim([0 0.5])
        ylabel('Probability Density','FontWeight','bold');
        title(sprintf('Cluster %d, %s-%d → %s-1', clmN, d_type{typeN,1}, scaleN, d_type{typeN,2}), ...
            'FontSize', 18, 'FontWeight', 'bold');

        saveas(fig, fullfile(pdDir, sprintf("Type%d_Cluster%d_PD.png", typeN, clmN)));
        close(fig)
    end
end
toc
