%% Setup & load data
close all; clear; clc;

requiredFile = "SouthAsia_DI.mat";
if ~isfile(requiredFile)
    error('Required file "%s" was not found in the current folder.', requiredFile);
end
if ~exist('select_best_copula_builtin', 'file')
    error('The helper function "select_best_copula_builtin.m" must be available on the MATLAB path.');
end

load(requiredFile)

requiredVars = ["SPI","SMI","SRI","gldas_lat","gldas_lon","gldas_land","shp_world"];
missingVars = requiredVars(~ismember(requiredVars, string(who)));
if ~isempty(missingVars)
    error('SouthAsia_DI.mat is missing required variables: %s', strjoin(cellstr(missingVars), ', '));
end

dataMap = struct('SPI', {SPI}, 'SMI', {SMI}, 'SRI', {SRI});

%% Output settings
resultFile = "DP_Grid.mat";
figureDir = fullfile("figures", "stationary");
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

%% Drought propagation model settings
margin_dist = 'Kernel';
thr_c = -1;
thr_e = -1;
copula_cand = {'Gaussian', 't', 'Clayton', 'Frank', 'Gumbel'};

%% Run analysis for each propagation pathway
tmp_di = {'SPI','SMI'; 'SPI','SRI'; 'SMI','SRI'};
tmp_header = {'SPI→SMI','SPI→SRI','SMI→SRI'};
tmp_scale = 1:6;

tic
PROB = cell(numel(tmp_header),1);
TS = cell(numel(tmp_header),1);

for i = 1:numel(tmp_header)
    PROB{i,1} = nan(size(gldas_land,1), size(gldas_land,2));
    TS{i,1} = nan(size(gldas_land,1), size(gldas_land,2));

    for latN = 1:size(gldas_lat,1)
        for lonN = 1:size(gldas_lon,2)
            fprintf('Analyzing %s [%d/%d, %d/%d]\n', tmp_header{i}, latN, size(gldas_lat,1), lonN, size(gldas_lon,2));

            if ~gldas_land(latN,lonN)
                continue
            end

            cause_matrix = nan(size(dataMap.(tmp_di{i,1}){1,1}, 3), numel(tmp_scale));
            for k = 1:numel(tmp_scale)
                cause_matrix(:,k) = squeeze(dataMap.(tmp_di{i,1}){k,1}(latN,lonN,:));
            end
            cause_matrix = cause_matrix(12:end,:);

            effect_vector = squeeze(dataMap.(tmp_di{i,2}){1,1}(latN,lonN,12:end));
            n_data = numel(effect_vector);

            if all(isnan(effect_vector)) || all(isnan(cause_matrix(:,1)))
                continue
            end

            try
                pd_effect = fitdist(effect_vector, margin_dist);
                v_effect = cdf(pd_effect, effect_vector);
            catch
                continue
            end

            results_family = cell(numel(tmp_scale), 1);
            results_params = cell(numel(tmp_scale), 1);
            results_aic = nan(numel(tmp_scale), 1);
            results_cond_prob = nan(numel(tmp_scale), 1);

            for ts = 1:numel(tmp_scale)
                current_cause = cause_matrix(:, ts);

                try
                    pd_cause = fitdist(current_cause, margin_dist);
                    u_cause = cdf(pd_cause, current_cause);

                    [family, params, aic] = select_best_copula_builtin([u_cause, v_effect], copula_cand);
                    results_family{ts} = family;
                    results_params{ts} = params;
                    results_aic(ts) = aic;

                    prob_cause_margin = sum(current_cause <= thr_c) / n_data;
                    prob_effect_margin = sum(effect_vector <= thr_e) / n_data;

                    if prob_cause_margin > 0
                        joint_prob = copulacdf(family, [prob_effect_margin, prob_cause_margin], params{:});
                        results_cond_prob(ts) = joint_prob / prob_cause_margin;
                    end
                catch
                    results_family{ts} = '';
                    results_params{ts} = {};
                    results_aic(ts) = nan;
                    results_cond_prob(ts) = nan;
                end
            end

            [max_prob, prop_duration] = max(results_cond_prob);
            if all(isnan(results_cond_prob))
                continue
            end

            thr_ts = 0.02;
            idx_df = find(abs(diff(results_cond_prob)) <= thr_ts, 1, 'first');
            if ~isempty(idx_df) && ~isnan(results_cond_prob(idx_df+1))
                max_prob = results_cond_prob(idx_df+1);
                prop_duration = idx_df + 1;
            end

            PROB{i,1}(latN,lonN) = max_prob;
            TS{i,1}(latN,lonN) = prop_duration;
        end
    end
end
toc

save(resultFile, "PROB", "TS", "tmp_di", "tmp_header", "tmp_scale", "thr_c", "thr_e", "margin_dist", "copula_cand", "-v7.3");

%% Plot map figure
tmp_al = string(('a':'z').').';
alN = 0;

tic
fig = figure('Units', 'pixels', 'Position', [0 0 850 1000], 'Color', 'w');
tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for headN = 1:numel(tmp_header)
    alN = alN + 1;
    nexttile; hold on;
    geoshow([shp_world.Y], [shp_world.X], 'DisplayType', 'polygon', 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'none');
    contourm(gldas_lat, gldas_lon, PROB{headN,1}, 'Fill', 'on', 'LineStyle', 'none');
    geoshow([shp_world.Y], [shp_world.X], 'DisplayType', 'line', 'LineWidth', 0.8, 'Color', [0.2 0.2 0.2]);
    xlim([60 98]); ylim([5.5 38.5]);
    clim([0 1]); colormap(gca, flipud(hot(10)));
    xlabel("Longitude",'FontWeight','bold'); ylabel("Latitude",'FontWeight','bold');
    if mod(alN,2) == 0, set(gca,'YLabel',[],'YTick',[]); end
    if ~ismember(alN,[5 6]), set(gca,'XLabel',[],'XTick',[]); end
    set(gca,"FontName","Arial","FontSize",12)
    title(sprintf("(%s) Propagation Probability %s", tmp_al(alN), tmp_header{headN}), 'FontSize', 14, 'FontWeight', 'bold')
    hold off

    alN = alN + 1;
    nexttile; hold on;
    geoshow([shp_world.Y], [shp_world.X], 'DisplayType', 'polygon', 'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'none');
    contourm(gldas_lat, gldas_lon, TS{headN,1}, 'Fill', 'on', 'LineStyle', 'none');
    geoshow([shp_world.Y], [shp_world.X], 'DisplayType', 'line', 'LineWidth', 0.8, 'Color', [0.2 0.2 0.2]);
    xlim([60 98]); ylim([5.5 38.5]);
    clim([0.5 6.5]);
    custom_map = [0.8 0.9 1; 0.6 0.8 1; 0.4 0.7 0.9; 1 0.8 0.4; 0.9 0.5 0.2; 0.8 0.2 0.1];
    colormap(gca, custom_map);
    xlabel("Longitude",'FontWeight','bold'); ylabel("Latitude",'FontWeight','bold');
    if mod(alN,2) == 0, set(gca,'YLabel',[],'YTick',[]); end
    if ~ismember(alN,[5 6]), set(gca,'XLabel',[],'XTick',[]); end
    set(gca,"FontName","Arial","FontSize",12)
    title(sprintf("(%s) Propagation Timescale %s", tmp_al(alN), tmp_header{headN}), 'FontSize', 14, 'FontWeight', 'bold')
    hold off
end
toc

saveas(fig, fullfile(figureDir, "DP_stationary_maps.png"));
