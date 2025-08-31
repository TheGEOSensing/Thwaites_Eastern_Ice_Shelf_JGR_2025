%% PINNING POINT ANALYSIS

% ----------------------- %
% Karen's data on MODIS
% ----------------------- % 
cd '/Users/student/Desktop/Objective 1/Clipped_Datasets_for_time_series_analysis'/'Karen''s Data'/

elon_files = dir('Elong_pp*.tif');
vv_files = dir('vv_pp*.tif');


for i = 1:size(vv_files)

% -------------------------------- %
% the vv files as velocity image
% -------------------------------- %
    image1 = vv_files(i).name;
    image = readgeoraster(image1); image(image == 0) = nan;
    vv_pp_MODIS(:,:,i) = double(image);

% -------------------------------- %
% elon files- longitudinal strain rates
% -------------------------------- %
    image1 = elon_files(i).name;
    image = readgeoraster(image1); image(image == 0) = nan; 
    elon_pp_MODIS(:,:,i) = double(image);
    
end

% -----------------
% number of elements
% ------------------
n_pp_MODIS = sum(~isnan(image(:)));
root_n_1_pp_MODIS = sqrt(n_pp_MODIS - 1);



%%
% ----------------------- %
% ITS_LIVE DATA
% ----------------------- % 
cd '/Users/student/Desktop/Objective 1/Clipped_Datasets_for_time_series_analysis/ITS_LIVE'

elon_files = dir('Elong_pp*.tif');
vv_files = dir('vv_pp*.tif');


for i = 1:size(vv_files)

% -------------------------------- %
% the vv files as velocity image
% -------------------------------- %
    image1 = vv_files(i).name;
    image = readgeoraster(image1); image(image == 0) = nan;
    vv_pp_ITS_LIVE(:,:,i) = double(image)/365;

% ------------------------------------- %
% elon files- longitudinal strain rates
% ------------------------------------- %
    image1 = elon_files(i).name;
    image = readgeoraster(image1); image(image == 0) = nan;
    elon_pp_ITS_LIVE(:,:,i) = double(image);
    
end

% -----------------
% number of elements
% ------------------
n_pp_ITS_LIVE = sum(~isnan(image(:)));
root_n_1_pp_ITS_LIVE = sqrt(n_pp_ITS_LIVE - 1);



clearvars -except  vv_pp_MODIS elon_pp_MODIS ...
    vv_pp_ITS_LIVE elon_pp_ITS_LIVE ...
    root_n_1_pp_ITS_LIVE root_n_1_pp_MODIS ...
    n_pp_ITS_LIVE n_pp_MODIS



%% Compute Time series of Elon and Velocity
% ------------------------------------------
pp_vv_time_series_MODIS = squeeze(nanmean(nanmean(vv_pp_MODIS,1),2));
pp_vv_time_series_ITS_LIVE = squeeze(nanmean(nanmean(vv_pp_ITS_LIVE,1),2));
total_pp_vv_time_series = [pp_vv_time_series_MODIS(2:13);pp_vv_time_series_ITS_LIVE(2:end)];


pp_elon_time_series_MODIS = squeeze(nansum(nansum(elon_pp_MODIS,1),2));
pp_elon_time_series_ITS_LIVE = squeeze(nansum(nansum(elon_pp_ITS_LIVE,1),2));
total_pp_elon_time_series = [pp_elon_time_series_MODIS(2:13); pp_elon_time_series_ITS_LIVE(2:end)];


% standard error calculation
% --------------------------
std_dev_pp_elon_MODIS = (squeeze(nanstd(elon_pp_MODIS, [], 1:2))) / root_n_1_pp_MODIS * n_pp_MODIS;
std_dev_pp_elon_ITS_LIVE = (squeeze(nanstd(elon_pp_ITS_LIVE, [], 1:2))) / root_n_1_pp_ITS_LIVE * n_pp_ITS_LIVE;
std_dev_pp_elon = [std_dev_pp_elon_MODIS(2:13); std_dev_pp_elon_ITS_LIVE(2:end)];


std_dev_pp_vv_MODIS = (squeeze(nanstd(vv_pp_MODIS, [], 1:2)));
std_dev_pp_vv_ITS_LIVE = (squeeze(nanstd(vv_pp_ITS_LIVE, [], 1:2)));
std_dev_pp_vv = [std_dev_pp_vv_MODIS(2:13); std_dev_pp_vv_ITS_LIVE(2:end)]/4;



%% CALCULATION OF COMPRESSION CHANGE AT THE PINNING POINT UPSTREAM AREA
% --------------------------------------------------------------------- %
years = 2002:2022;
figure(3)

% --------------
% velocity 
% --------------
yyaxis left
errorbar(years,total_pp_vv_time_series,std_dev_pp_vv,'LineStyle','-.','LineWidth',2,'color',"k");
ylabel('Ice Flow Speed (m day^-^1)','fontsize',16,'fontweight','bold'); 
ylim([0 1.45])

% --------------------------------------------------------------------- %
% let's make a patch of the years for identifying different time period
% --------------------------------------------------------------------- %
patch([2001 2007 2007 2001], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','facealpha',0.15,'edgecolor','none')
patch([2007 2012 2012 2007], [max(ylim) max(ylim) min(ylim) min(ylim)],'g','facealpha',0.15,'edgecolor','none')
patch([2012 2017 2017 2012], [max(ylim) max(ylim) min(ylim) min(ylim)],'y','facealpha',0.15,'edgecolor','none')
patch([2017 2023 2023 2017], [max(ylim) max(ylim) min(ylim) min(ylim)],'r','facealpha',0.15,'edgecolor','none')

% --------------------------
% longitudinal strain rates
% --------------------------
yyaxis right

% --------------------------------------------------------------------- %
% make the bars red for the negative (compressive) strain and blue for
% positive (extension) strain rates
% --------------------------------------------------------------------- %

hold on
for i = 1:length(total_pp_elon_time_series)
    h=bar(years(i),total_pp_elon_time_series(i));

    if total_pp_elon_time_series(i) < 0
        set(h,'FaceColor','r','facealpha',0.7);
    elseif total_pp_elon_time_series(i) > 0
        set(h,'FaceColor','#0072BD','facealpha',0.7);
    end
end

er = errorbar(years,total_pp_elon_time_series,...
    std_dev_pp_elon ,'linewidth',2,'color','#7E2F8E');    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylim([-0.025  0.025])
ylabel('Total Longitudinal strain rates (day^-^1)','fontsize',16,'fontweight','bold');

hold on 

% ------------------------------------------ %
% labelling the x-axis of the graph
% ------------------------------------------ %
xticks(2002:2022)
xlim([2001.5 2022.5])

% ------------------------------------------ %
% design the axis of the plot
% ------------------------------------------ %
ax = gca; ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;

ax.YAxis(2).Color = "r";
ax.YAxis(1).Color = 'k';

set(gca, 'SortMethod', 'depth')

