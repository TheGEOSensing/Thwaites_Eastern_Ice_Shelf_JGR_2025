cd('/Users/rishi/Desktop/Objective 1/Thwaites_fracture_Sentinel/Length_TEIS_fractures')


% Use a for loop to get te length of the fractures from digitized shapefile
files = dir('Length*.shp');

for i = 1:size(files)
    file = files(i).name;
    fracture = shaperead(file);
    fracture_length = nansum([fracture(:).Length_lin]);
    fracture_lengths(i) = fracture_length; 
end

err_fracture_lengths = fracture_lengths * 0.0329; % calculated from Sentinel-1 fractures


%% Use a for loop to get the area of the internal melanges from digitized shapefile

cd '/Users/rishi/Desktop/Objective 1/Thwaites_fracture_Sentinel/Thwaites_Sentinel_Internal_melange/TEIS_Melange_areas'

files = dir('Area*.shp');

for i = 1:size(files)
    file = files(i).name;
    melange = shaperead(file);
    area = sum([melange(:).Area]);
    melange_area(i) = area; 
end

err_melange_area = melange_area* 0.0957; % calculated from landsat fractures



%% Create a datetime array from January 2020 to December 2022

monthly_dates = datetime(2020, 1, 1):calmonths(1): datetime(2022, 12, 31);;

%%


% ------------------------------------------
% Left axis for the shear fracture length
% ------------------------------------------
yyaxis left
plot(monthly_dates, fracture_lengths,'LineWidth',2.5,'Color','k')
ylabel('Length of Fracture (km)','fontsize',16,'fontweight','bold')
hold on
errorbar(monthly_dates, fracture_lengths, err_fracture_lengths,'k','LineStyle','none','LineWidth',1.25);
ylim([160 235])

% ------------------------------------------
% Right axis for the internal melange areas
% ------------------------------------------
yyaxis right
b = bar(monthly_dates, melange_area);
set(b,'FaceColor','r','facealpha',0.6);
ylabel('Area of Internal Melange (sq. km.)','fontsize',16,'fontweight','bold')
hold on 
errorbar(monthly_dates, melange_area, err_melange_area,'k','LineStyle','none','LineWidth',1.25);
xlim([monthly_dates(1) monthly_dates(31)]);
ylim([1 16])

patch([monthly_dates(7) monthly_dates(10) monthly_dates(10) monthly_dates(7)], [max(ylim) max(ylim) min(ylim) min(ylim)],'k','facealpha',0.15,'edgecolor','none')


% ---------------------------
% Design axis for the plots 
% ---------------------------
ax = gca; ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;
ax.YAxis(2).Color = "r";
ax.YAxis(1).Color = 'k';


%%

% ------------------------------------------
% Left axis for the shear fracture length
% ------------------------------------------
yyaxis left
hold on

% Convert datetime to datenum
xnum = datenum(monthly_dates);

% Compute bounds
upper = fracture_lengths + err_fracture_lengths;
lower = fracture_lengths - err_fracture_lengths;

% Prepare fill region
x_fill = [xnum, fliplr(xnum)];
y_fill = [upper, fliplr(lower)];

% Plot filled envelope
fill(x_fill, y_fill, 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.13);

% Plot main curve (also using datenum)
plot(xnum, fracture_lengths, 'k', 'LineWidth', 2.5);

% Fix x-axis formatting
datetick('x', 'mmm-yyyy', 'keeplimits');
ylabel('Length of Fracture (km)', 'fontsize', 16, 'fontweight', 'bold');
ylim([160 235])


% ------------------------------------------
% Right axis for the internal melange areas
% ------------------------------------------
yyaxis right
hold on

% Convert datetime to datenum for x-axis consistency
xnum = datenum(monthly_dates);

% Plot bar chart
b = bar(xnum, melange_area);
set(b, 'FaceColor', 'r', 'FaceAlpha', 0.6);
errorbar(xnum, melange_area, err_melange_area, 'k', ...
         'LineStyle', 'none', 'LineWidth', 1.25);
ylabel('Area of Internal Melange (sq. km.)', ...
       'fontsize', 16, 'fontweight', 'bold');
ylim([1 16])

% Patch region for a date range
patch([xnum(7) xnum(10) xnum(10) xnum(7)], ...
      [max(ylim) max(ylim) min(ylim) min(ylim)], ...
      'k', 'facealpha', 0.1, 'edgecolor', 'none');

% Set x-axis limits and formatting
xlim([xnum(1) xnum(31)]);
datetick('x', 'mmm-yyyy', 'keeplimits');

% ---------------------------
% Design axis for the plots 
% ---------------------------
ax = gca; ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;
ax.YAxis(2).Color = "r";
ax.YAxis(1).Color = 'k';
