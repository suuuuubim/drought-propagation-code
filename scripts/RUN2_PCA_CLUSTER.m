%% Setup & load data
close all; clear; clc;

requiredFile = "SouthAsia_DI.mat";
if ~isfile(requiredFile)
    error('Required file "%s" was not found in the current folder.', requiredFile);
end

load(requiredFile)

requiredVars = ["SPI","SMI","SRI","gldas_land"];
missingVars = requiredVars(~ismember(requiredVars, string(who)));
if ~isempty(missingVars)
    error('SouthAsia_DI.mat is missing required variables: %s', strjoin(cellstr(missingVars), ', '));
end

figureDir = fullfile("figures", "clustering");
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

%% Cluster grids based on drought-index time series
nTime = size(SPI{1},3) - 12;
nScale = 6;
gridN = size(gldas_land,1) * size(gldas_land,2);
tmp_x = nan(gridN, nTime * 3 * nScale);

idx = 1;
for i = 1:size(SPI{1},1)
    for k = 1:size(SPI{1},2)
        tmp_feature = [];
        flag_skip = false;
        for scN = 1:nScale
            tmp_spi = squeeze(SPI{scN}(i,k,13:end));
            tmp_smi = squeeze(SMI{scN}(i,k,13:end));
            tmp_sri = squeeze(SRI{scN}(i,k,13:end));

            if any(isnan(tmp_spi)) || any(isnan(tmp_smi)) || any(isnan(tmp_sri))
                flag_skip = true;
                break;
            end
            tmp_feature = [tmp_feature; tmp_spi; tmp_smi; tmp_sri];
        end

        if ~flag_skip
            tmp_x(idx,:) = tmp_feature';
        end
        idx = idx + 1;
    end
end

tmp_valid = all(~isnan(tmp_x),2);
tmp_xvalid = tmp_x(tmp_valid,:);

%% PCA
[~, tmp_score, ~, ~, tmp_expl] = pca(tmp_xvalid);
tmp_cumExpl = cumsum(tmp_expl);
tmp_nComp = find(tmp_cumExpl >= 90, 1);
tmp_xnew = tmp_score(:,1:tmp_nComp);

%% K-means clustering
rng(1,'twister');
clusterN = 8;
[idx_c, tmp_c] = kmeans(tmp_xnew, clusterN, 'Distance', 'correlation', 'Replicates', 5);

cluster_map = nan(size(SPI{1},1), size(SPI{1},2));
idx_grid = 1;
for i = 1:size(SPI{1},1)
    for k = 1:size(SPI{1},2)
        if tmp_valid(idx_grid)
            cluster_map(i,k) = idx_c(sum(tmp_valid(1:idx_grid)));
        end
        idx_grid = idx_grid + 1;
    end
end

%% Plot cluster map
tmp_color = summer(max(cluster_map(:)));
[~, idx_tmp] = sort(rand(max(cluster_map(:)),1));
tmp_color = tmp_color(idx_tmp,:);

fig = figure('Color','w'); hold on;
cluster_label = cluster_map;
cluster_label(isnan(cluster_label)) = 0;
stats = regionprops(cluster_label, 'Centroid');

imagesc(cluster_map, 'AlphaData', ~isnan(cluster_map))
colormap(tmp_color);

for cN = 1:max(cluster_map(:))
    mask = cluster_label == cN;
    if any(mask(:))
        centroid = mean([stats(cN).Centroid], 1);
        text(centroid(1), centroid(2), num2str(cN), "FontSize", 15, ...
            "FontWeight", "bold", "Color", "k", "FontName", "Times New Roman", ...
            "HorizontalAlignment", "center");
    end
end

title('Cluster Map')
set(gca,'YDir','reverse','FontName','Times New Roman','FontSize',12)
ax = gca;
ax.YAxis.Visible = 'off';
ax.XAxis.Visible = 'off';
hold off

save("NS_Result_C8.mat", "cluster_map", "clusterN", "idx_c", "tmp_c", "tmp_valid", "tmp_nComp", "tmp_cumExpl", "-v7.3");
saveas(fig, fullfile(figureDir, "cluster_map_C8.png"));
