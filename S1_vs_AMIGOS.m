cd '/Users/rishi/Desktop/Objective 1/Adrian''s data on TEIS monthly velocity'
load Adrian_S1_monthly_data+AMIGOS.mat


%% -------------------------------------------------------------------------------- %
% plot the time series of the shear strain rates for the pinning point shear zone
% -------------------------------------------------------------------------------- %
time = datetime(2020, 1, 1):calmonths(1): datetime(2022, 11, 1);

figure
yyaxis left
plot(time, shear_mean_eshear_time_series_2019_2021,'LineWidth',2.5,'LineStyle','-.')
hold on
errorbar(time, shear_mean_eshear_time_series_2019_2021, std_err_shear_boundary_eshear_2013_2022,...
    'linewidth',2,'linestyle','none','Color','k')

% manage the plot distributions and style
% ----------------------------------------
ax = gca; 
ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;
ylabel('Average Shear Strain Rates (day^-^1)','fontsize',14,'fontweight','bold');
xlim([time(1), datetime(2022, 7, 1)]);


% Plot the trenline by calculating coefficients
% ----------------------------------------------
coeffs = polyfit((1:35),shear_mean_eshear_time_series_2019_2021,1);
shear_fit = polyval(coeffs,(1:35));

hold on
plot(time, shear_fit,'linestyle','--', 'LineWidth',1.5,'color','k')


% put the r-squared value in the plot from the model
% --------------------------------------------------
mdl = fitlm((1:35),shear_boundary_eshear_time_series_2019_2021,1);
text(datetime(2021, 03, 01), 1.13*10^-6,['Adjusted R² = ', sprintf('%.2f', mdl.Rsquared.Adjusted)],...
    'fontweight','bold' ,'Color', '#0072BD', 'FontSize', 12,'HorizontalAlignment', 'left')

text(datetime(2021, 03, 01), 1.19*10^-6, ['Shearing Trend = ',sprintf('%.2f x10^-^7 day^-^2', coeffs(1)*10^7)],...
    'fontweight','bold' ,'Color', '#0072BD', 'FontSize', 12);
ylim([0.1*10^-6  1.25*10^-6])



% yy axis for the divergence over the shear zone
% --------------------------------------------------
yyaxis right

plot(time,shear_boundary_div_time_series_2019_2021,'LineWidth',2.5,'LineStyle','-.')
hold on
errorbar(time, shear_boundary_div_time_series_2019_2021, std_err_shear_boundary_div_2013_2022,...
    'linewidth',2,'linestyle','none','Color','k')
ylabel('Average Flow Divergence (day^-^1)','fontsize',14,'fontweight','bold');
ylim([0 7*10^-6])


% Plot the trenline by calculating coefficients
% ----------------------------------------------
coeffs1 = polyfit((1:35),shear_boundary_div_time_series_2019_2021,1);
div_fit = polyval(coeffs1,(1:35));

hold on
plot(time, div_fit,'linestyle','--', 'LineWidth',1.5,'color','k')


mdl1 = fitlm((1:35),shear_boundary_div_time_series_2019_2021,1);
text(datetime(2021, 03, 01), 0.3*10^-6,['Adjusted R² = ', sprintf('%.2f', mdl1.Rsquared.Adjusted)],...
    'fontweight','bold' ,'Color', '#D95319', 'FontSize', 12,'HorizontalAlignment', 'left')

text(datetime(2021, 03, 01), 0.65*10^-6, ['Divergence Trend = ',sprintf('%.2f x10^-^7 day^-^2', coeffs1(1)*10^6)],...
    'fontweight','bold' ,'Color', '#D95319', 'FontSize', 12);


patch([DOY_2020_may(1) DOY_2020_may(end) DOY_2020_may(end) DOY_2020_may(1)]...
    , [max(ylim) max(ylim) min(ylim) min(ylim)], 'k','facealpha',0.15,'edgecolor','none')

patch([DOY_2021_feb(1) DOY_2022_march(end) DOY_2022_march(end) DOY_2021_feb(1)]...
    , [max(ylim) max(ylim) min(ylim) min(ylim)], 'k','facealpha',0.15,'edgecolor','none')



%% AMIGOS GPS data plotting over shear strain rates %%%
% ------------------------------------------------------

figure
yyaxis right
scatter(DOY2020,Speedmd,'.','r'); 
hold on
plot(DOY2020, Speedmd, 'r','LineStyle','--'); 



% quarters of the time
% -------------------------
% plot(DOY_2020_may, AMIGOS_2020_may, 'r','LineStyle','-.'); 
scatter(DOY_2020_may, AMIGOS_2020_may,'.','r');
text(datetime(2020, 06, 20), 1.75,[sprintf('%.2f mm/day^2', slope1 )],... % linear trend
    'fontweight','bold' ,'Color', 'r', 'FontSize', 14);


 scatter(DOY_2021_july, AMIGOS_2021_july,'.','r');
 text(datetime(2020, 11, 01), 1.86,[sprintf('%.2f mm/day^2', slope2 )],... % linear trend
     'fontweight','bold' ,'Color', 'r', 'FontSize', 14);
xlim([time(1), datetime(2022, 7, 1)]);


% design the axis and the patches
% -------------------------------
ax = gca; 
ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;
ax.YAxis(2).Color = "r";

ylabel('AMIGOS GPS Ice Flow Speed (m/day)','fontsize',14,'fontweight','bold')
patch([DOY_2020_may(1) DOY_2020_may(end) DOY_2020_may(end) DOY_2020_may(1)]...
    , [max(ylim) max(ylim) min(ylim) min(ylim)], 'k','facealpha',0.15,'edgecolor','none')

patch([DOY_2021_feb(1) DOY_2022_march(end) DOY_2022_march(end) DOY_2021_feb(1)]...
    , [max(ylim) max(ylim) min(ylim) min(ylim)], 'k','facealpha',0.15,'edgecolor','none')



% plot the time series of the shear strain rates for the pinning point shear zone
% -------------------------------------------------------------------------------- %
time = datetime(2020, 1, 1):calmonths(1): datetime(2022, 11, 1);


yyaxis left
plot(time, shear_mean_eshear_time_series_2019_2021,'LineWidth',2.5,'LineStyle','-.')
hold on
errorbar(time, shear_mean_eshear_time_series_2019_2021, std_err_shear_boundary_eshear_2013_2022,...
    'linewidth',2,'linestyle','none','Color','k')


% manage the plot distributions and style
% ----------------------------------------
ax = gca; 
ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;
ylabel('Average Shear Strain Rates (day^-^1)','fontsize',14,'fontweight','bold');
xlim([time(1), datetime(2022, 7, 1)]);


% Plot the trenline by calculating coefficients
% ----------------------------------------------
coeffs = polyfit((1:35),shear_mean_eshear_time_series_2019_2021,1);
shear_fit = polyval(coeffs,(1:35));

hold on
plot(time, shear_fit,'linestyle','--', 'LineWidth',1.5,'color','k')

% put the r-squared value in the plot from the model
% --------------------------------------------------
mdl = fitlm((1:35),shear_boundary_eshear_time_series_2019_2021,1);
text(datetime(2021, 03, 01), 1.13*10^-6,['Adjusted R² = ', sprintf('%.2f', mdl.Rsquared.Adjusted)],...
    'fontweight','bold' ,'Color', '#0072BD', 'FontSize', 12,'HorizontalAlignment', 'left')

text(datetime(2021, 03, 01), 1.19*10^-6, ['Shearing Trend = ',sprintf('%.2f x10^-^7 day^-^2', coeffs(1)*10^7)],...
    'fontweight','bold' ,'Color', '#0072BD', 'FontSize', 12);
ylim([0.1*10^-6  1.25*10^-6])