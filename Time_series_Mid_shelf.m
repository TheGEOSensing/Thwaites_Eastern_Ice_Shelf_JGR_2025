
%% MID_SHELF ANALYSIS

% ----------------------- %
% Karen's data on MODIS
% ----------------------- % 
cd '/Users/student/Desktop/Objective 1/Clipped_Datasets_for_time_series_analysis'/'Karen''s Data'/

elon_files = dir('Elong_mid*.tif');
vv_files = dir('vv_clip*.tif');


for i = 1:size(vv_files)

% -------------------------------- %
% the vv files as velocity image
% -------------------------------- %
    image1 = vv_files(i).name;
    image = readgeoraster(image1); image(image == 0) = nan;
    vv_mid_shelf_MODIS(:,:,i) = double(image);

% -------------------------------- %
% elon files- longitudinal strain rates
% -------------------------------- %
    image1 = elon_files(i).name;
    image = readgeoraster(image1); image(image == 0) = nan; 
    elon_mid_shelf_MODIS(:,:,i) = double(image);

end

% -----------------
% number of elements
% ------------------
n_mid_shelf_MODIS = sum(~isnan(image(:)));
root_n_1_mid_shelf_MODIS = sqrt(n_mid_shelf_MODIS - 1);


clearvars -except  vv_mid_shelf_MODIS elon_mid_shelf_MODIS ...
    vv_mid_shelf_ITS_LIVE elon_mid_shelf_ITS_LIVE ...
    root_n_1_mid_shelf_ITS_LIVE root_n_1_mid_shelf_MODIS ...
    n_mid_shelf_ITS_LIVE n_mid_shelf_MODIS 

%%
% ----------------------- %
% ITS_LIVE DATA
% ----------------------- % 
cd '/Users/student/Desktop/Objective 1/Clipped_Datasets_for_time_series_analysis/ITS_LIVE'

elon_files = dir('Elong_mid*.tif');
vv_files = dir('vv_mid*.tif');


for i = 1:size(vv_files)

    % -------------------------------- %
    % the vv files as velocity image
    % -------------------------------- %
    image1 = vv_files(i).name;
    image = readgeoraster(image1); image(image == 0) = nan;
    vv_mid_shelf_ITS_LIVE(:,:,i) = double(image)/365;

    % ------------------------------------- %
    % elon files- longitudinal strain rates
    % ------------------------------------- %
    image1 = elon_files(i).name;
    image = readgeoraster(image1); image(image == 0) = nan;
    elon_mid_shelf_ITS_LIVE(:,:,i) = double(image);
    
end

% -----------------
% number of elements
% ------------------
n_mid_shelf_ITS_LIVE = sum(~isnan(image(:)));
root_n_1_mid_shelf_ITS_LIVE = sqrt(n_mid_shelf_ITS_LIVE - 1);



clearvars -except  vv_mid_shelf_MODIS elon_mid_shelf_MODIS ...
    vv_mid_shelf_ITS_LIVE elon_mid_shelf_ITS_LIVE ...
    root_n_1_mid_shelf_ITS_LIVE root_n_1_mid_shelf_MODIS ...
    n_mid_shelf_ITS_LIVE n_mid_shelf_MODIS



%% Compute Time series of Elon and Velocity

mid_shelf_vv_time_series_MODIS = squeeze(nanmean(nanmean(vv_mid_shelf_MODIS,1),2));
mid_shelf_vv_time_series_ITS_LIVE = squeeze(nanmean(nanmean(vv_mid_shelf_ITS_LIVE,1),2));
total_mid_shelf_vv_time_series = [mid_shelf_vv_time_series_MODIS(1:13);mid_shelf_vv_time_series_ITS_LIVE(3:end)];


mid_shelf_elon_time_series_MODIS = squeeze(nansum(nansum(elon_mid_shelf_MODIS,1),2));
mid_shelf_elon_time_series_ITS_LIVE = squeeze(nansum(nansum(elon_mid_shelf_ITS_LIVE,1),2));
total_mid_shelf_elon_time_series = [mid_shelf_elon_time_series_MODIS(1:13);mid_shelf_elon_time_series_ITS_LIVE(3:end)];


% standard error calculation
% --------------------------
std_err_mid_shelf_elon_MODIS = squeeze(nanstd(elon_mid_shelf_MODIS, [], 1:2)/root_n_1_mid_shelf_MODIS) * n_mid_shelf_MODIS;
std_err_mid_shelf_elon_ITS_LIVE = squeeze(nanstd(elon_mid_shelf_ITS_LIVE, [], 1:2)/root_n_1_mid_shelf_ITS_LIVE) * n_mid_shelf_ITS_LIVE;
std_err_mid_shelf_elon = [std_err_mid_shelf_elon_MODIS(1:13); std_err_mid_shelf_elon_ITS_LIVE(3:end)];


std_err_mid_shelf_vv_MODIS = squeeze(nanstd(vv_mid_shelf_MODIS, [], 1:2) / root_n_1_mid_shelf_MODIS);
std_err_mid_shelf_vv_ITS_LIVE = squeeze(nanstd(vv_mid_shelf_ITS_LIVE, [], 1:2) / root_n_1_mid_shelf_ITS_LIVE);
std_err_mid_shelf_vv = [std_err_mid_shelf_vv_MODIS(1:13); std_err_mid_shelf_vv_ITS_LIVE(3:end)];


%%

years = 2002:2022;
figure(2)

% --------------
% velcoity 
% --------------
yyaxis left
plot(years,total_mid_shelf_vv_time_series,'LineStyle','-.','LineWidth',3,'color','k');
hold on
a = errorbar(years,total_mid_shelf_vv_time_series,std_err_mid_shelf_vv,'LineStyle','none',...
    'LineWidth',2,'color','k');
ylabel('Ice flow Speed (m/day)','fontsize',16,'fontweight','bold');
ylim([0 3.5])


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
hold on
for i = 1:length(total_mid_shelf_elon_time_series)
    h = bar(years(i),total_mid_shelf_elon_time_series(i));

    if total_mid_shelf_elon_time_series(i) < 0
        set(h,'FaceColor','r','facealpha',0.7);
    elseif total_mid_shelf_elon_time_series(i) > 0
        set(h,'FaceColor',"#0072BD",'facealpha',0.8);
    end
end


% ------------------------------------------ %
% include the errorbar in the bar diagram
% ------------------------------------------ %
er = errorbar(years,total_mid_shelf_elon_time_series,std_err_mid_shelf_elon,'linewidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

ylabel('Longitudinal strain rates (day^-^1)','fontsize',16,'fontweight','bold');
ylim([-6.2*10^-3  6*10^-3])


% ------------------------------------------ %
% labelling the x-axis of the graph
% ------------------------------------------ %
xticks(2002:2022)
xlim([2001.4 2022.6])


% ------------------------------------------ %
% design the axis of the plot
% ------------------------------------------ %
ax = gca; ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;

ax.YAxis(2).Color = "r";
ax.YAxis(1).Color = 'k';

set(gca, 'SortMethod', 'depth')
