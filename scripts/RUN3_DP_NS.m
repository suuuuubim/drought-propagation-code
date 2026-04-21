%% Setup & load data
close all; clear; clc;

requiredFiles = ["SouthAsia_DI.mat", "NS_Result_C8.mat"];
for f = requiredFiles
    if ~isfile(f)
        error('Required file "%s" was not found in the current folder.', f);
    end
end

load("SouthAsia_DI.mat")
load("NS_Result_C8.mat")

requiredVars = ["SPI","SMI","SRI","cluster_map"];
missingVars = requiredVars(~ismember(requiredVars, string(who)));
if ~isempty(missingVars)
    error('Required variables are missing: %s', strjoin(cellstr(missingVars), ', '));
end

dataMap = struct('SPI', {SPI}, 'SMI', {SMI}, 'SRI', {SRI});

outputDir = fullfile("figures", "DroughtPropagation_NS");
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% Drought propagation model settings
thr_c = -1;
thr_e = -1;
options = optimoptions('fminunc','Display','off');
pad_len = 24;

%% Copula selection and fitting per climate zone
d_type = {'SPI','SMI'; 'SPI','SRI'; 'SMI','SRI'};

tic
cond_prob_results = cell(size(d_type,1),1);

for typeN = 1:size(d_type,1)
    tic
    for clmN = 1:max(cluster_map(:))
        [latN, lonN] = find(cluster_map == clmN);
        cond_prob_all = [];
        cond_prob_ts = [];

        for grdN = 1:numel(latN)
            cause_matrix = [];
            for i = 1:6
                tmp_cause = squeeze(dataMap.(d_type{typeN,1}){i,1}(latN(grdN),lonN(grdN),13:end));
                cause_matrix = [cause_matrix tmp_cause];
            end
            effect_vector = squeeze(dataMap.(d_type{typeN,2}){1,1}(latN(grdN),lonN(grdN),13:end));

            if sum(isnan(cause_matrix),'all') == 0 && sum(isnan(effect_vector)) == 0
                n_scales = size(cause_matrix,2);
                n_obs = size(cause_matrix,1);

                cause_padded = [flipud(cause_matrix(1:pad_len,:)); cause_matrix; flipud(cause_matrix(end-pad_len+1:end,:))];
                effect_padded = [flipud(effect_vector(1:pad_len)); effect_vector; flipud(effect_vector(end-pad_len+1:end))];
                time_padded = (1:length(cause_padded))'/12;

                cond_prob_ns = nan(n_obs, n_scales);

                for tsN = 1:n_scales
                    fprintf("Cluster %d: (%d/%d)th grid, Non-stationary fitting for %s-%d\n", ...
                        clmN, grdN, numel(latN), d_type{typeN,1}, tsN);

                    current_cause = cause_padded(:,tsN);

                    k = 4;
                    knots = linspace(min(time_padded), max(time_padded), floor(max(time_padded)/10) + 2);
                    t_aug = [repmat(knots(1),1,k-1), knots, repmat(knots(end),1,k-1)];
                    B_spline = spcol(t_aug, k, time_padded);
                    n_basis = size(B_spline,2);

                    negLogL_effect = @(p) -sum(log(normpdf(effect_padded, B_spline*p(1:n_basis), exp(B_spline*p(n_basis+1:end)))));
                    p_effect = fminunc(negLogL_effect, [zeros(n_basis,1); log(std(effect_padded))*ones(n_basis,1)], options);
                    mu_e_t = B_spline*p_effect(1:n_basis);
                    sigma_e_t = exp(B_spline*p_effect(n_basis+1:end));

                    negLogL_cause = @(p) -sum(log(normpdf(current_cause, B_spline*p(1:n_basis), exp(B_spline*p(n_basis+1:end)))));
                    p_cause = fminunc(negLogL_cause, [zeros(n_basis,1); log(std(current_cause))*ones(n_basis,1)], options);
                    mu_c_t = B_spline*p_cause(1:n_basis);
                    sigma_c_t = exp(B_spline*p_cause(n_basis+1:end));

                    u_ns = normcdf(current_cause, mu_c_t, sigma_c_t);
                    v_ns = normcdf(effect_padded, mu_e_t, sigma_e_t);
                    u_ns = min(max(u_ns,1e-4),1-1e-4);
                    v_ns = min(max(v_ns,1e-4),1-1e-4);

                    negLogL_copula = @(p) -sum(arrayfun(@(t) log(copulapdf('Gaussian', [u_ns(t), v_ns(t)], ...
                        min(max(tanh(B_spline(t,:)*p), -0.999), 0.999)) + 1e-12), 1:length(u_ns)));
                    p_copula = fminunc(negLogL_copula, zeros(n_basis,1), options);
                    rho_t = tanh(B_spline*p_copula);
                    rho_t = min(max(rho_t, -0.999), 0.999);

                    mu_c_t = mu_c_t(pad_len+1:end-pad_len);
                    sigma_c_t = sigma_c_t(pad_len+1:end-pad_len);
                    mu_e_t = mu_e_t(pad_len+1:end-pad_len);
                    sigma_e_t = sigma_e_t(pad_len+1:end-pad_len);
                    rho_t = rho_t(pad_len+1:end-pad_len);

                    u_thr_t = normcdf(thr_c, mu_c_t, sigma_c_t);
                    v_thr_t = normcdf(thr_e, mu_e_t, sigma_e_t);

                    for obsN = 1:n_obs
                        if u_thr_t(obsN) > 1e-6
                            joint_prob = mvncdf([norminv(u_thr_t(obsN)), norminv(v_thr_t(obsN))], [0 0], [1, rho_t(obsN); rho_t(obsN), 1]);
                            cond_prob_ns(obsN, tsN) = joint_prob / u_thr_t(obsN);
                        end
                    end
                end

                mean_prob = mean(cond_prob_ns, 'omitmissing');
                [tmp_max, ~] = max(mean_prob);
                idx_new = find(abs(mean_prob - tmp_max) <= 0.01);

                if ~isempty(idx_new)
                    cond_prob_all = [cond_prob_all cond_prob_ns(:,idx_new(1))];
                    cond_prob_ts = [cond_prob_ts; idx_new(1)];
                end
            end
        end

        cond_prob_results{typeN,1}{clmN,1} = cond_prob_all;
        cond_prob_results{typeN,1}{clmN,2} = cond_prob_ts;

        if isempty(cond_prob_all)
            warning('Cluster %d for type %d returned no valid results.', clmN, typeN);
            continue
        end

        n_obs_plot = size(cond_prob_all, 1);

        fig = figure('Position',[200 200 1000 300], 'Color', 'w'); hold on;
        tmp_qt = quantile(cond_prob_all',0.10,1);
        tmp_qt = [tmp_qt; quantile(cond_prob_all',0.25,1)];
        tmp_qt = [tmp_qt; quantile(cond_prob_all',0.50,1)];
        tmp_qt = [tmp_qt; quantile(cond_prob_all',0.75,1)];
        tmp_qt = [tmp_qt; quantile(cond_prob_all',0.90,1)];
        fill([1:n_obs_plot fliplr(1:n_obs_plot)], [tmp_qt(2,:) fliplr(tmp_qt(4,:))], [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'Edgecolor', 'none')
        fill([1:n_obs_plot fliplr(1:n_obs_plot)], [tmp_qt(1,:) fliplr(tmp_qt(5,:))], [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'Edgecolor', 'none')
        plot(tmp_qt(3,:), 'Color', [1 1 1]*0.2, 'LineWidth', 2);

        tmp_p = polyfit(1:n_obs_plot, tmp_qt(3,:), 1);
        tmp_f = polyval(tmp_p, 1:n_obs_plot);
        plot(1:n_obs_plot, tmp_f, 'b--', 'LineWidth', 2)

        text(max(1, n_obs_plot-120), 0.1, sprintf("%0.4fx+%0.2f", tmp_p(1), tmp_p(2)), ...
            'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 15, 'Color', 'b')
        grid on; grid minor; box on;
        xticks([1:12*5:min(600,n_obs_plot) n_obs_plot]); xlim([1 n_obs_plot]); ylim([0 1]);
        xticklabels(string([1976 1980:5:2020 2024]))
        xlabel("Year"); ylabel(sprintf("P(%s|%s)", d_type{typeN,2}, d_type{typeN,1}));
        legend("Quantile 25-75%", "Quantile 10-90%", "Median", "Linear Fit", "location", "north", "Orientation", "horizontal")
        title(sprintf("Drought Propagation Probability in Cluster %d", clmN))
        set(gca,"FontName","Times New Roman","FontWeight","bold","FontSize",12)
        hold off;

        saveas(fig, fullfile(outputDir, sprintf("Type%d_Cluster%d_quantile.png", typeN, clmN)));
    end
    toc
end
toc

save("cond_prob_results.mat", "cond_prob_results", "-v7.3");

%% Plot timescale histogram per cluster
alphabet = string(('a':'z').').';

tic
for typeN = 1:size(d_type,1)
    fig = figure('Position', [100 50+250*(typeN-1) 1400 400], 'Color', 'w');
    tiledlayout(2, 4, 'Padding', 'tight', 'TileSpacing', 'compact');

    for clmN = 1:max(cluster_map(:))
        nexttile; hold on;
        if isempty(cond_prob_results{typeN,1}{clmN,2})
            text(0.5, 0.5, 'No valid data', 'HorizontalAlignment', 'center');
            axis off
            continue
        end

        hs = histogram(cond_prob_results{typeN,1}{clmN,2});
        hs_c = hs.Values;
        hs_e = hs.BinEdges;
        [~, hs_idx] = max(hs_c);
        histogram('BinEdges', hs_e(hs_idx:hs_idx+1), 'BinCounts', hs_c(hs_idx), 'FaceColor', 'r')
        xticks(1:1:6); xlim([0.2 6.8]);
        set(gca,'FontName','Arial','FontSize',14)
        if clmN > 4
            xlabel("Timescale (month)",'FontWeight','bold');
        end
        if ismember(clmN,[1 5])
            ylabel("Frequency",'FontWeight','bold');
        end
        title(sprintf("(%s) Cluster %d", alphabet(clmN), clmN), 'FontWeight', 'bold', 'FontSize', 20)
        box off; hold off;
    end

    saveas(fig, fullfile(outputDir, sprintf("Type%d_cluster_timescale_hist.png", typeN)));
end
toc

%% Plot quantile series for all clusters
tic
for typeN = 1:size(d_type,1)
    fig = figure('Position', [0+1000*(typeN-1) 0 1200 1000], 'Color', 'w');
    tiledlayout(4, 2, 'Padding', 'tight', 'TileSpacing', 'loose');

    for clmN = 1:max(cluster_map(:))
        nexttile; hold on;

        tmp_data = cond_prob_results{typeN,1}{clmN,1};
        if isempty(tmp_data)
            text(0.5, 0.5, 'No valid data', 'HorizontalAlignment', 'center');
            axis off
            continue
        end

        n_obs_plot = size(tmp_data,1);
        tmp_qt = quantile(tmp_data',0.10,1);
        tmp_qt = [tmp_qt; quantile(tmp_data',0.25,1)];
        tmp_qt = [tmp_qt; quantile(tmp_data',0.50,1)];
        tmp_qt = [tmp_qt; quantile(tmp_data',0.75,1)];
        tmp_qt = [tmp_qt; quantile(tmp_data',0.90,1)];
        fill([1:n_obs_plot fliplr(1:n_obs_plot)], [tmp_qt(2,:) fliplr(tmp_qt(4,:))], [0.5 0.5 0.5], 'FaceAlpha', 0.4, 'Edgecolor', 'none')
        fill([1:n_obs_plot fliplr(1:n_obs_plot)], [tmp_qt(1,:) fliplr(tmp_qt(5,:))], [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'Edgecolor', 'none')
        plot(tmp_qt(3,:), 'Color', [1 1 1]*0.2, 'LineWidth', 2);

        tmp_p = polyfit(1:n_obs_plot, tmp_qt(3,:), 1);
        tmp_f = polyval(tmp_p, 1:n_obs_plot);
        plot(1:n_obs_plot, tmp_f, 'b--', 'LineWidth', 2)
        text(max(1, n_obs_plot-140), 0.1, sprintf("%0.4fx+%0.2f", tmp_p(1), tmp_p(2)), ...
            'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 15, 'Color', 'b')

        xticks([1:12*5:min(600,n_obs_plot) n_obs_plot]); xlim([1 n_obs_plot]); ylim([0 1]);
        xticklabels(string([1976 1980:5:2020 2024]))
        if ismember(clmN,[7 8])
            xlabel("Year","FontWeight","bold");
        end

        mode_scale = mode(cond_prob_results{typeN,1}{clmN,2}, 'all');
        ylabel(sprintf("P(%s-1|%s-%d)", d_type{typeN,2}, d_type{typeN,1}, mode_scale), "FontWeight", "bold");
        set(gca,"FontName","Arial","FontSize",12)
        title(sprintf("(%s) Drought Propagation in Cluster %d", alphabet(clmN), clmN), "FontWeight", "bold", "FontSize", 18)
        hold off;
    end

    saveas(fig, fullfile(outputDir, sprintf("Type%d_all_clusters_quantile.png", typeN)));
end
toc

%% Mann-Kendall test
tbl_type = [];
tbl_cluster = [];
tbl_senSlope = [];
tbl_Z = [];
tbl_p = [];
tbl_sig = string([]);

tic
mk_results = cell(size(d_type,1),1);

for typeN = 1:size(d_type,1)
    for clmN = 1:max(cluster_map(:))
        tmp_data = cond_prob_results{typeN,1}{clmN,1};
        if isempty(tmp_data)
            continue
        end

        tmp_qt = quantile(tmp_data', 0.5, 1);

        tmp_s = 0;
        tmp_cnt = 0;
        tmp_slope = nan(length(tmp_qt)*(length(tmp_qt)-1)/2,1);

        for i = 1:length(tmp_qt)-1
            for j = i+1:length(tmp_qt)
                tmp_s = tmp_s + sign(tmp_qt(j)-tmp_qt(i));
                tmp_cnt = tmp_cnt + 1;
                tmp_slope(tmp_cnt) = (tmp_qt(j)-tmp_qt(i))/(j-i);
            end
        end

        tmp_uniqueX = unique(tmp_qt);
        if length(tmp_qt) == length(tmp_uniqueX)
            tmp_vars = length(tmp_qt)*(length(tmp_qt)-1)*(2*length(tmp_qt)+5)/18;
        else
            tmp_tp = zeros(length(tmp_uniqueX),1);
            for k = 1:length(tmp_uniqueX)
                tmp_tp(k) = sum(tmp_qt == tmp_uniqueX(k));
            end
            tmp_vars = (length(tmp_qt)*(length(tmp_qt)-1)*(2*length(tmp_qt)+5) - sum(tmp_tp.*(tmp_tp-1).*(2*tmp_tp+5))) / 18;
        end

        if tmp_s > 0
            tmp_z = (tmp_s-1)/sqrt(tmp_vars);
        elseif tmp_s < 0
            tmp_z = (tmp_s+1)/sqrt(tmp_vars);
        else
            tmp_z = 0;
        end

        tmp_pval = 2*(1-normcdf(abs(tmp_z),0,1));
        tmp_h = abs(tmp_z) > norminv(1-0.05/2,0,1);
        tmp_senSlope = median(tmp_slope,'omitmissing');

        mk_results{typeN,1}{clmN,1}.h = tmp_h;
        mk_results{typeN,1}{clmN,1}.p = tmp_pval;
        mk_results{typeN,1}{clmN,1}.Z = tmp_z;
        mk_results{typeN,1}{clmN,1}.S = tmp_s;
        mk_results{typeN,1}{clmN,1}.VarS = tmp_vars;
        mk_results{typeN,1}{clmN,1}.sen_slope = tmp_senSlope;
        mk_results{typeN,1}{clmN,1}.median_ts = tmp_qt;

        tbl_type(end+1,1) = typeN;
        tbl_cluster(end+1,1) = clmN;
        tbl_senSlope(end+1,1) = tmp_senSlope;
        tbl_Z(end+1,1) = tmp_z;
        tbl_p(end+1,1) = tmp_pval;

        if tmp_pval < 0.001
            tbl_sig(end+1,1) = "***";
        elseif tmp_pval < 0.01
            tbl_sig(end+1,1) = "**";
        elseif tmp_pval < 0.05
            tbl_sig(end+1,1) = "*";
        else
            tbl_sig(end+1,1) = "ns";
        end
    end
end
toc

MK_table = table(tbl_type, tbl_cluster, tbl_senSlope, tbl_Z, tbl_p, tbl_sig, ...
    'VariableNames', {'Type','Cluster','SenSlope','Z','p_value','Significance'});

save("MK_results.mat", "mk_results", "MK_table", "-v7.3");
writetable(MK_table, "MK_results.csv");
