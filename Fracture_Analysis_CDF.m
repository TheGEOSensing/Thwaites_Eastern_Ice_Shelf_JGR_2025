
cd '/Users/rishi/Desktop/Objective 1/Clipped_fractures_TEIS/Fracture_Orientation_epsg3031'


% Read a year shapefile and get the fracture orientation from that
fracture_2002 = shaperead('fractures_2002.shp');
Length_2002 = [fracture_2002(:).length];

fracture_2006 = shaperead('fractures_2006.shp');
Length_2006 = [fracture_2006(:).length];

fracture_2010 = shaperead('fractures_2010.shp');
Length_2010 = [fracture_2010(:).length];

fracture_2014 = shaperead('fractures_2013.shp');
Length_2014 = [fracture_2014(:).length];

fracture_2018 = shaperead('fractures_2018.shp');
Length_2018 = [fracture_2018(:).length];

fracture_2022 = shaperead('fractures_2021.shp');
Length_2022 = [fracture_2022(:).length];

% Take a threshold for the lower limit (David's Suggestion)
% -------------------------------------------------------------------
Length_2002 = Length_2002(Length_2002 >= 0.5 );
Length_2006 = Length_2006(Length_2006 >= 0.5 );
Length_2010 = Length_2010(Length_2010 >= 0.5 );
Length_2014 = Length_2014(Length_2014 >= 0.5 );
Length_2018 = Length_2018(Length_2018 >= 0.5 );
Length_2022 = Length_2022(Length_2022 >= 0.5 );


% combined histograom operation for different years
% ---------------------------------------------------
max_value = max([Length_2002, Length_2006, Length_2010, Length_2014, Length_2018, Length_2022]);
edges = [0:2:8, max_value];

% prepare the histogram counts to plot
% -------------------------------------
h1 = histcounts(Length_2002,edges);
h2 = histcounts(Length_2006,edges);
h3 = histcounts(Length_2010,edges);
h4 = histcounts(Length_2014,edges);
h5 = histcounts(Length_2018,edges);
h6 = histcounts(Length_2022,edges);


%%


figure
clf

% yyaxis left
yyaxis left
% ------------------------------------------------------------ %
% plot the bar graphs with the fracture frequency distribution %
% ------------------------------------------------------------ %
b = bar(edges(1:end-1), [h1; h2; h3; h4; h5; h6]');

% Apply colors to bars
colors = (brewermap(6,"RdBu"));
for i = 1:6
    b(i).FaceColor = colors(i,:);
    b(i).EdgeColor = 'k';
    b(i).LineWidth = 1.3;
end

% --------------------
% Design the legend
% --------------------
l = legend(' 2002', ' 2006', ' 2010', ' 2014', ' 2018', ' 2021');
l.FontSize = 14; l.EdgeColor = 'w';

% -----------------------------------
% design the attributes of the figure
% -----------------------------------
ax1 = gca; ax1.FontSize = 14;
ax1.FontWeight = 'bold';
ax1.LineWidth = 2.5; grid off; box on;

ylabel('Number of Fractures');
xlim([-1 9]); ylim([0 70]);
xticklabels(ax1,{'\leq 2 km','2.01 - 4 km','4.01 - 6 km', '6.01 - 8 km', '> 8 km'})

% -------------------------------------------------------------- %
% prepare yy axis for the cumulative distribution of the histogram
% -------------------------------------------------------------- %

yyaxis right
cdfplot(Length_2002);
hold on;
cdfplot(Length_2006);
cdfplot(Length_2010);
cdfplot(Length_2014);
cdfplot(Length_2018);
cdfplot(Length_2022);

% -----------------------------------
% Set different colors for each line
% ------------------------------------
colors = flipud(brewermap(6, 'RdBu'));
lines = findobj(gca, 'Type', 'Line');

for i = 1:numel(lines)
    set(lines(i), 'Color', colors(i,:), 'LineWidth', 2,'LineStyle','-','Marker','none');
end

% -----------------------
% Make the figure labels
% -----------------------
xlabel('Fracture Length Distributions');
ylabel('Cumulative Frequency');

legend('2002', '2006', '2010', '2014', '2018', '2021','Location','east');
hold off;
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = 'k';
grid off

% Create secondary x-axis on top
ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', 'YAxisLocation', 'right', 'Color', 'none');

% Set the x-ticks and x-tick labels for the top axis
ax2.XTick = ax1.XTick;
ax2.XTickLabel = ax1.XTickLabel;

% Make the top x-axis labels visible
ax2.XTickLabel = { '0 km','2 km','4 km','6 km', '8 km'};
ax2.YTickLabel = [] ;

% Set the x limits to be the same
ax2.XLim = ax1.XLim;
set(ax2, 'color', 'none', 'box', 'off');
linkaxes([ax1 ax2], 'x');
ax2.TickLength = [0 0];
xlabel(ax2, 'Fracture length for CDF');
ax2.FontWeight = 'bold'; 
ax2.FontSize = 14;
ax2.LineWidth = 2.5; grid off; box on;

