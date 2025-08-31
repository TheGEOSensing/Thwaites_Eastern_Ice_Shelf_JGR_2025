%% PINNING POINT ANALYSIS

% ----------------------- %
% Karen's data on MODIS
% ----------------------- % 
cd '/Users/student/Desktop/Objective 1/Clipped_Datasets_for_time_series_analysis'/'Karen''s Data'/

div_files = dir('div_shear*.tif');
eshear_files = dir('Eshear*.tif');


for i = 1:size(eshear_files)

% -------------------------------- %
% eshear files over shear zone
% -------------------------------- %
    image1 = eshear_files(i).name;
    image = abs(readgeoraster(image1)); 
    image(image == 0) = nan; % outside the shear zone
    eshear_MODIS(:,:,i) = double(image);
    n_sz_MODIS_eshear = sum(~isnan(image(:)));

% -------------------------------- %
% divergence files- divergence
% -------------------------------- %
    image1 = div_files(i).name;
    image = readgeoraster(image1); 
    image(image <= 0) = nan;    % only the positive values (divergence) and outside
    div_MODIS(:,:,i) = double(image);
    n_sz_MODIS_div = sum(~isnan(image(:)));
end

% -----------------
% number of elements
% ------------------
root_n_1_sz_MODIS_eshear = sqrt(n_sz_MODIS_eshear - 1);
root_n_1_sz_MODIS_div = sqrt(n_sz_MODIS_div - 1);


%%
% ----------------------- %
% ITS_LIVE DATA
% ----------------------- % 
cd '/Users/student/Desktop/Objective 1/Clipped_Datasets_for_time_series_analysis/ITS_LIVE'

div_files = dir('div*.tif');
eshear_files = dir('Eshear*.tif');


for i = 1:size(eshear_files)

% -------------------------------- %
% eshear files over shear zone
% -------------------------------- %
    image1 = eshear_files(i).name;
    image = abs(readgeoraster(image1)); 
    image(image == 0) = nan;
    eshear_ITS_LIVE(:,:,i) = double(image);
    n_sz_ITS_LIVE_eshear = sum(~isnan(image(:)));

% -------------------------------- %
% divergence files- divergence
% -------------------------------- %
    image1 = div_files(i).name;
    image = readgeoraster(image1); 
    image(image <= 0) = nan; 
    div_ITS_LIVE(:,:,i) = double(image);
    n_sz_ITS_LIVE_div = sum(~isnan(image(:)));

end

% -----------------
% number of elements
% ------------------
root_n_1_sz_ITS_LIVE_eshear = sqrt(n_sz_ITS_LIVE_eshear - 1);
root_n_1_sz_ITS_LIVE_div = sqrt(n_sz_ITS_LIVE_div - 1);


% ------------------------------------------
%% Compute Time series of Elon and Velocity
% ------------------------------------------
sz_eshear_time_series_MODIS = squeeze(nansum(nansum(eshear_MODIS,1),2));
sz_eshear_time_series_ITS_LIVE = squeeze(nansum(nansum(eshear_ITS_LIVE,1),2));
total_sz_eshear_time_series = [sz_eshear_time_series_MODIS(2:12); sz_eshear_time_series_ITS_LIVE(1:end)];


sz_div_time_series_MODIS = squeeze(nansum(nansum((div_MODIS),1),2));
sz_div_time_series_ITS_LIVE = squeeze(nansum(nansum((div_ITS_LIVE),1),2));
total_sz_div_time_series = [sz_div_time_series_MODIS(2:12); sz_div_time_series_ITS_LIVE(1:end)];


% --------------------------
% standard error calculation
% --------------------------
std_err_sz_eshear_MODIS = (squeeze(nanstd(eshear_MODIS, [], 1:2)) / root_n_1_sz_MODIS_eshear) * n_sz_MODIS_eshear;
std_err_sz_eshear_ITS_LIVE = (squeeze(nanstd(eshear_ITS_LIVE, [], 1:2)) / root_n_1_sz_ITS_LIVE_eshear) * n_sz_ITS_LIVE_eshear;
std_err_sz_eshear = [std_err_sz_eshear_MODIS(2:12); std_err_sz_eshear_ITS_LIVE(1:end)];


std_err_sz_div_MODIS = (squeeze(nanstd(div_MODIS, [], 1:2)) / root_n_1_sz_MODIS_div) * n_sz_MODIS_div;
std_err_sz_div_ITS_LIVE = (squeeze(nanstd(div_ITS_LIVE, [], 1:2)) / root_n_1_sz_ITS_LIVE_div) * n_sz_ITS_LIVE_div;
std_err_sz_div = [std_err_sz_div_MODIS(2:12); std_err_sz_div_ITS_LIVE(1:end)];



%% Plot the shear zone figure
% ----------------------------
years = 2002:2022;
figure(4)


yyaxis left
errorbar(years, total_sz_eshear_time_series, std_err_sz_eshear,...
    'Linestyle','-.', 'LineWidth',2.5,'color','k');
ylabel('Total Shear Strain Rate (day^-^1)','fontsize',16,'fontweight','bold');
hold on

% --------------------------------------------------------------------- %
% let's make a patch of the years for identifying different time period
% --------------------------------------------------------------------- %
patch([2001 2007 2007 2001], [max(ylim) max(ylim) min(ylim) min(ylim)],'c','facealpha',0.15,'edgecolor','none')
patch([2007 2012 2012 2007], [max(ylim) max(ylim) min(ylim) min(ylim)],'g','facealpha',0.15,'edgecolor','none')
patch([2012 2017 2017 2012], [max(ylim) max(ylim) min(ylim) min(ylim)],'y','facealpha',0.15,'edgecolor','none')
patch([2017 2023 2023 2017], [max(ylim) max(ylim) min(ylim) min(ylim)],'r','facealpha',0.15,'edgecolor','none')


yyaxis right
errorbar(years, total_sz_div_time_series, std_err_sz_div,...
    'Linestyle','-.', 'LineWidth',2,'color','#0072BD');
ylabel('Total Ice Flow Divergence (day^-^1)','fontsize',16,'fontweight','bold');


% ----------------------------------- %
% labelling the x-axis of the graph
% ----------------------------------- %
xticks(2002:2022)
xlim([2001 2023])

% ----------------------------- %
% design the axis of the plot
% ----------------------------- %
ax = gca; ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = '#0072BD';


