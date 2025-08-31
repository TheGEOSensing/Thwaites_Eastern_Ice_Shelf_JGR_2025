cd '/Users/rishi/Desktop/Objective 1/Adrian''s data on TEIS monthly velocity'
load Adrian_S1_monthly_data.mat


% Data Processed and saved in the .mat file, Don't need to run the strain rates again 

%% get all the velocity products and calculate the strain rates from that
% ----------------------------------------------------------------------- %
vv_files = dir('*mag.tif');
vx_files = dir('*x.tif');
vy_files = dir('*y.tif');

for i = 1:size(vv_files)

% -------------------------------- %
% the vv files as velocity image
% -------------------------------- %
    image1 = vv_files(i).name;
    image = readgeoraster(image1);
    vv(:,:,i) = double(image);

% -------------------------------- %
% vx files as the longitudinal velocity
% -------------------------------- %
    image1 = vx_files(i).name;
    image = readgeoraster(image1);
    vx(:,:,i) = double(image);
    
% -------------------------------- %
% vy files as the transverse velocity
% -------------------------------- %
    image1 = vy_files(i).name;
    [image,R] = readgeoraster(image1);
    vy(:,:,i) = double(image);
end


% ------------------------------------------------------- %
% calculate the strain rates from the velocity components
% ------------------------------------------------------- %
for i = 1:size(vx,3)
    [elon,etrans,eshear,eEff] = ...
        strainRates( vx(:,:,i), vy(:,:,i), 500);
    % ------------------------------------------------------- %
    % capture strain rates in individual 3d matrices
    % ------------------------------------------------------- %
    eps_long(:,:,i) = elon;
    eps_trans(:,:,i) = etrans;
    eps_shear(:,:,i) = eshear;
    eps_Eff(:,:,i) = eEff;    
end

%% get the latitude and longitude from the image file
% ------------------------------------------------------- %
% longitude extent 
xmin = R.XWorldLimits(1);  
xmax = R.XWorldLimits(2);

% latitiude extent
ymin = R.YWorldLimits(1);  
ymax = R.YWorldLimits(2);  

% get the X and Y from the max and min of longitude (divide it with the
% size of the raster
Y = ymin: (ymax - ymin)/size(vv,1): ymax;
X = xmin: (xmax - xmin)/size(vv,2): xmax;

Y = imresize(Y,[1,size(vv,1)]);
X = imresize(X,[1,size(vv,2)]);

% meshgrid - latitude and longitude and then convert it to lat and long
[X,Y] = meshgrid(X, Y);
[Latitude, Longitude] = ps2ll(X,Y); % polar stereographic to lat long

%%
% 
% figure
% pcolorps(Latitude,Longitude, flipud(abs(eps_shear(:,:,2)))); caxis([-0.00005 0.00005])
% hold on
% colormap(brewermap([],"RdBu"));
% 
% %%
% figure
% pcolorps(Latitude,Longitude, flipud(abs(vv(:,:,71))));
% colormap(brewermap([],"Blues"));

%%
cd('/Users/student/Desktop/Objective 1/Thwaites_Rectangle/Shear_Boundary_shapefiles/')
shapefile = shaperead('Shea_zone_smaller.shp'); % read shapefile

% ------------------------------------------------------- %
% Get the spatial information from the shapefile
% ------------------------------------------------------- %
ry = shapefile.Y(1:end-1);
rx = shapefile.X(1:end-1);

% ------------------------------------------------------- %
% clip or mask the raster data from shapefile
% this is a seperate script file (inpsquad.m)
% ------------------------------------------------------- %
mask = inpolygon(Latitude,Longitude,ry,rx);


for i = 1:size(vx,3)-1
    % --------------------------------------------------- %
    % make a masked for each of the datasets
    % multiply the masked with the matrix to clip
    % --------------------------------------------------- %
    masked_ehsear = flipud(eps_shear(:,:,i)).* mask;      

    % --------------------------------------------------- %
    % combine the temporal datasets together
    % --------------------------------------------------- %
    shear_boundary_eshear(:,:,i) = abs(masked_ehsear);
end

% shearing time series over shear zone %
% --------------------------------------
shear_boundary_eshear_time_series_2019_2021 = squeeze(nansum(nansum(abs(shear_boundary_eshear),1),2));
shear_mean_eshear_time_series_2019_2021 = squeeze(nanmean(nanmean(abs(shear_boundary_eshear),1),2));


% number of sample calculation %
% -------------------------------
n_shear_zone = (size(shear_boundary_eshear,1) * size(shear_boundary_eshear,2));

% square root of n-1 for errorbars %
% ----------------------------------
root_n_1_shear_zone = sqrt((size(shear_boundary_eshear,1) * size(shear_boundary_eshear,2)) - 1);

% calculate the error value from the standard devaition using 
% ------------------------------------------------------------ %
std_err_shear_boundary_eshear_2013_2022 = squeeze(nanstd(abs(shear_boundary_eshear), [], 1:2)/root_n_1_shear_zone)*n_shear_zone;
 

%% 
% -------------------------------------------------------------------------------- %
% plot the time series of the shear strain rates for the pinning point shear zone
% -------------------------------------------------------------------------------- %
time = datetime(2020, 1, 1):calmonths(1): datetime(2022, 11, 1);

figure
yyaxis left
plot(time,shear_boundary_eshear_time_series_2019_2021,'LineWidth',2.5,'LineStyle','-.')
hold on
errorbar(time, shear_boundary_eshear_time_series_2019_2021, std_err_shear_boundary_eshear_2013_2022,...
    'linewidth',1.5,'linestyle','none','Color','k')

% manage the plot distributions and style
% ----------------------------------------
ax = gca; 
ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.LineWidth = 2.5; grid off; box on;
ylabel('Total shear strain rate (day^-^1)','fontsize',14,'fontweight','bold');
xlim([time(1), datetime(2022, 8, 1)]);
ylim([0.01 0.24])

% Plot the trenline by calculating coefficients
% ----------------------------------------------
coeffs = polyfit((1:35),shear_boundary_eshear_time_series_2019_2021,1);
shear_fit = polyval(coeffs,(1:35));

hold on
plot(time, shear_fit,'linestyle','-.', 'LineWidth',1.5,'color','k')


% put the r-squared value in the plot from the model
% --------------------------------------------------
mdl = fitlm((1:35),shear_boundary_eshear_time_series_2019_2021,1);
text(datetime(2021, 05, 1), 0.06,['Adjusted R² = ', sprintf('%.2f', mdl.Rsquared.Adjusted)],...
    'fontweight','bold' ,'Color', '#0072BD', 'FontSize', 14)

text(datetime(2021, 05, 1), 0.072, ['Shearing Trend = ',sprintf('%.3f day^-^2', mdl.Coefficients.Estimate)],...
    'fontweight','bold' ,'Color', '#0072BD', 'FontSize', 14);


%%% AMIGOS GPS data plotting over shear strain rates %%%
% --------------------------------------------------------
yyaxis right
scatter(DOY2020,Speedmd,'.','r'); 
hold on
plot(DOY2020, Speedmd, 'r','LineStyle','--'); 

% quarters of the time
% -------------------------
% plot(DOY_2020_may, AMIGOS_2020_may, 'r','LineStyle','-.'); 
scatter(DOY_2020_may, AMIGOS_2020_may,'.','r');
text(datetime(2020, 06, 20), 1.75,[sprintf('%.2f mm/day^2', 0.36 )],... % linear trend
    'fontweight','bold' ,'Color', 'r', 'FontSize', 14);


 scatter(DOY_2021_july, AMIGOS_2021_july,'.','r');
 text(datetime(2021, 03, 01), 2.3,[sprintf('%.2f mm/day^2', 0.92 )],... % linear trend
     'fontweight','bold' ,'Color', 'r', 'FontSize', 14);


% design the axis and the patches
% -------------------------------
ylabel('AMIGOS GPS Ice Flow Speed (m/day)','fontsize',14,'fontweight','bold')
ax.YAxis(2).Color = "r";
patch([DOY_2020_may(1) DOY_2020_may(end) DOY_2020_may(end) DOY_2020_may(1)]...
    , [max(ylim) max(ylim) min(ylim) min(ylim)], 'k','facealpha',0.15,'edgecolor','none')

patch([DOY_2021_feb(1) DOY_2022_march(end) DOY_2022_march(end) DOY_2021_feb(1)]...
    , [max(ylim) max(ylim) min(ylim) min(ylim)], 'k','facealpha',0.15,'edgecolor','none')




%%

%%% AMIGOS GPS data sgementation %%%
% -------------------------------- %

AMIGOS_2020_may = Speedmd (1 : 9575);
AMIGOS_2021_feb = Speedmd (9776 : 21407);
AMIGOS_2021_july = Speedmd (21408 : 33669);
AMIGOS_2022_march = Speedmd (37749 : end);


DOY_2020_may = DOY2020 (1 : 9575);
DOY_2021_feb = DOY2020 (9776 : 21407);
DOY_2021_july = DOY2020 (21408 : 33669);
DOY_2022_march = DOY2020 (37749 : end);


coeffs_2020_may = polyfit((1:9575), AMIGOS_2020_may,1);
coeffs_2021_feb = polyfit((9776 : 21407), AMIGOS_2021_feb,1);
coeffs_2021_july = polyfit((21408 : 33669), AMIGOS_2021_july,1);
coeffs_2022_march = polyfit((37749 : 49280), AMIGOS_2022_march,1);


polyfit_2020_may = polyval(coeffs_2020_may,(1 : 9575));
polyfit_2021_feb = polyval(coeffs_2021_feb,(9776 : 21407));
polyfit_2021_july = polyval(coeffs_2021_july,(21408 : 33669));
polyfit_2022_march = polyval(coeffs_2022_march,(37749 : 49280));


mdl_1 = fitlm((1:9575),AMIGOS_2020_may,1);
mdl_2 = fitlm((9776 : 21407),AMIGOS_2021_feb,1);
mdl_3 = fitlm((21408 : 33669),AMIGOS_2021_july,1);
mdl_4 = fitlm((37749 : 49280),AMIGOS_2022_march,1);


%% 