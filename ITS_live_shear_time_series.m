cd('/Users/student/Desktop/Objective 1/Its_live vleocity_image_pairs')
load('Thwaites_velocity_image_Pair_products.mat','velocity','vel_x','vel_y','Latitude','Longitude')

%%
for i = 1:size(vel_x,3)
  % -------------------------------------------------------- %
  % ------- Calculate STRAIN RATES for each of the year ---- %
  % -------------------------------------------------------- %
    % Karen's code on strain rate computation %

    [elon,etrans,eshear,exGrid,eyGrid,exyGrid,eEff] = ...
        strainRates( vel_x(:,:,i)/365, vel_y(:,:,i)/365, 500);

    % capture strain rates in individual 3d matrices
    eps_long(:,:,i) = elon;
    eps_trans(:,:,i) = etrans;
    eps_shear(:,:,i) = eshear;

    eps_xx(:,:,i) = exGrid;
    eps_yy(:,:,i) = eyGrid;
    eps_xy(:,:,i) = exyGrid;

    eps_Eff(:,:,i) = eEff;
end


%% Shapefile clipping operation
cd('/Users/student/Desktop/Objective 1/Thwaites_Rectangle/Shear_Boundary_shapefiles/')

shapefile = shaperead('Shear_boundary_rectangle.shp'); % read shapefile

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


for i = 1:size(vel_x,3)
    % ------------------------------------------------------------- %
    % make a masked for each of the datasets
    % multiply the masked with the matrix to clip
    % finally make the zero values as NaN values to get proper mean
    % ------------------------------------------------------------- %
    masked_vv = flipud(velocity(:,:,i)).* mask;     masked_vv(masked_vv == 0) = NaN;
    masked_vx = flipud(vel_x(:,:,i)).* mask;        masked_vx(masked_vx == 0) = NaN;
    masked_vy = flipud(vel_y(:,:,i)).* mask;        masked_vy(masked_vy == 0) = NaN;

    masked_elon = flipud(eps_long(:,:,i)).* mask;       masked_elon(masked_elon == 0) = NaN;
    masked_ehsear = flipud(eps_shear(:,:,i)).* mask;    masked_ehsear(masked_ehsear == 0) = NaN;
    masked_etrans = flipud(eps_trans(:,:,i)).* mask;    masked_etrans(masked_etrans == 0) = NaN;

    % --------------------------------------------------------- %
    % combine the temporal datasets together and make 3D matrix
    % --------------------------------------------------------- %
    shear_boundary_vv(:,:,i) = masked_vv;
    shear_boundary_vx(:,:,i) = masked_vx;
    shear_boundary_vy(:,:,i) = masked_vy;

    shear_boundary_elon(:,:,i) = masked_elon;
    shear_boundary_eshear(:,:,i) = masked_ehsear;
    shear_boundary_etrans(:,:,i) = masked_etrans;
end

n_shear_zone = sum(~isnan(masked_ehsear(:)));

% 
% figure
% pcolorps(Latitude,Longitude, (masked_ehsear(:,:,1)));
% 
% figure
% pcolorps(Latitude,Longitude, flipud(abs(eps_shear(:,:,1))));


%%
% --------------------------------------------------------------- %

shapefile = shaperead('Pinning_Point.shp'); % read shapefile

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


for i = 1:size(vel_x,3)
    % ------------------------------------------------------------- %
    % make a masked for each of the datasets
    % multiply the masked with the matrix to clip
    % finally make the zero values as NaN values to get proper mean
    % ------------------------------------------------------------- %
    masked_vv = flipud(velocity(:,:,i)).* mask;     masked_vv(masked_vv == 0) = NaN;
    masked_vx = flipud(vel_x(:,:,i)).* mask;        masked_vx(masked_vx == 0) = NaN;
    masked_vy = flipud(vel_y(:,:,i)).* mask;        masked_vy(masked_vy == 0) = NaN;

    masked_elon = flipud(eps_long(:,:,i)).* mask;       masked_elon(masked_elon == 0) = NaN;
    masked_ehsear = flipud(eps_shear(:,:,i)).* mask;    masked_ehsear(masked_ehsear == 0) = NaN;
    masked_etrans = flipud(eps_trans(:,:,i)).* mask;    masked_etrans(masked_etrans == 0) = NaN;


    % --------------------------------------------------- %
    % combine the temporal datasets together
    % --------------------------------------------------- %
    pinning_point_vv(:,:,i) = masked_vv;

    pinning_point_elon(:,:,i) = masked_elon;
    pinning_point_eshear(:,:,i) = masked_ehsear;
    pinning_point_etrans(:,:,i) = masked_etrans;
end

% number of elements
n_pinning_point = sum(~isnan(masked_elon(:)));

% 
% figure
% pcolorps(Latitude,Longitude, (masked_ehsear(:,:,1)));
% 
% figure
% pcolorps(Latitude,Longitude, flipud(abs(eps_shear(:,:,1))));


%% Mid Shelf dynamics calculation
% ------------------------------------------------------- %

shapefile = shaperead('Mid_shelf_Rectangle.shp'); % read shapefile

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


for i = 1:size(vel_x,3)
    % --------------------------------------------------- %
    % make a masked for each of the datasets
    % multiply the masked with the matrix to clip
    % --------------------------------------------------- %
    masked_vv = flipud(velocity(:,:,i)).* mask;     masked_vv(masked_vv == 0) = NaN;
    masked_vx = flipud(vel_x(:,:,i)).* mask;        masked_vx(masked_vx == 0) = NaN;
    masked_vy = flipud(vel_y(:,:,i)).* mask;        masked_vy(masked_vy == 0) = NaN;

    masked_elon = flipud(eps_long(:,:,i)).* mask;       masked_elon(masked_elon == 0) = NaN;
    masked_ehsear = flipud(eps_shear(:,:,i)).* mask;    masked_ehsear(masked_ehsear == 0) = NaN;
    masked_etrans = flipud(eps_trans(:,:,i)).* mask;    masked_etrans(masked_etrans == 0) = NaN;

    % --------------------------------------------------- %
    % combine the temporal datasets together
    % --------------------------------------------------- %
    mid_shelf_vv(:,:,i) = masked_vv;
    mid_shelf_vx(:,:,i) = masked_vx;
    mid_shelf_vy(:,:,i) = masked_vy;

    mid_shelf_elon(:,:,i) = masked_elon;
    mid_shelf_eshear(:,:,i) = masked_ehsear;
    mid_shelf_etrans(:,:,i) = masked_etrans;
end

n_mid_shelf = sum(~isnan(masked_ehsear(:)));

% figure
% pcolorps(Latitude,Longitude, (masked_ehsear(:,:,1)));
% 
% figure
% pcolorps(Latitude,Longitude, flipud(abs(eps_shear(:,:,1))));

%% ------------------------------------------------------- %
% shear zone time series of strain rates
% ------------------------------------------------------- %
shear_boundary_eshear_time_series_2013_2022 = squeeze(nansum(nansum(abs(shear_boundary_eshear),1),2));

% ------------------------------------------------------- %
% Mid shelf time series of velocity
% ------------------------------------------------------- %
mid_shelf_vv_time_series_2013_2022 = squeeze(nanmean(nanmean(mid_shelf_vv/365,1),2));


% ------------------------------------------------------- %
% Mid shelf time series of strain rates
% ------------------------------------------------------- %
mid_shelf_elon_time_series_2013_2022  = squeeze(nansum(nansum(mid_shelf_elon,1),2));


% ------------------------------------------------------- %
% Pinning Point time series of velocity
% ------------------------------------------------------- %
pinning_point_vv_time_series_2013_2022 = squeeze(nanmean(nanmean(pinning_point_vv/365,1),2));


% ------------------------------------------------------- %
% Pinning Point time series of strain rates
% ------------------------------------------------------- %
pinning_point_elon_time_series_2013_2022 = squeeze(nansum(nansum(pinning_point_elon,1),2));



%% ------------------------------------------------------- %
% Standard Deviation for the calculation of the errorbar 
% -------------------------------------------------------- %

% number of sample calculation %
% -------------------------------
% n_shear_zone = sum(~isnan(shear_boundary_vv(:)))/size(shear_boundary_vv,3);
% n_mid_shelf = sum(~isnan(mid_shelf_vv(:))) / size(mid_shelf_vv, 3);
% n_pinning_point = sum(~isnan(pinning_point_vv(:))) / size(pinning_point_vv, 3);

% ----------------------------------
% square root of n-1 for errorbars %
% ----------------------------------
root_n_1_shear_zone = sqrt(n_shear_zone - 1);
root_n_1_mid_shelf = sqrt(n_mid_shelf - 1);
root_n_1_pinning_point = sqrt(n_pinning_point - 1);

% ------------------------------------------------------------ %
% calculate the error value from the standard devaition using 
% ------------------------------------------------------------ %

% mid-shelf area
% --------------

std_err_mid_shelf_elon_2013_2022 = squeeze(nanstd(mid_shelf_elon, [], [1, 2])/root_n_1_mid_shelf)*n_mid_shelf; % because its sum
std_err_mid_shelf_vv_2013_2022 = squeeze(nanstd(mid_shelf_vv/365, [], [1, 2])/root_n_1_mid_shelf);

% shear boundary
% --------------
std_err_shear_boundary_eshear_2013_2022 = squeeze(nanstd(abs(shear_boundary_eshear), [], [1, 2])/root_n_1_shear_zone)*n_shear_zone;
 
% pinning point
% --------------
std_err_pinning_point_vv_2013_2022 = squeeze(nanstd(pinning_point_vv/365, [], [1, 2])/root_n_1_pinning_point); % beacuse it's a mean
std_err_pinning_point_elon_2013_2022 = squeeze(nanstd(pinning_point_elon, [], [1, 2])/root_n_1_pinning_point)*n_pinning_point;


% ------------------------------------------------------- %
%% Remove the unnecessary variables
% ------------------------------------------------------- %

clearvars -except ...
mid_shelf_vv_time_series_2000_2018...
mid_shelf_elon_time_series_2000_2018...
shear_boundary_eshear_time_series_2000_2018...
mid_shelf_vv_time_series_2013_2022 ...
mid_shelf_elon_time_series_2013_2022...
shear_boundary_eshear_time_series_2013_2022...
pinning_point_vv_time_series_2013_2022 ...
pinning_point_elon_time_series_2013_2022 ...
pinning_point_vv_time_series_2000_2018 ...
pinning_point_elon_time_series_2000_2018 ...
std_err_mid_shelf_elon_2013_2022 ...
std_err_mid_shelf_vv_2013_2022...
std_err_shear_boundary_eshear_2013_2022...
std_err_pinning_poing_elon_2013_2022...
std_err_pinning_point_vv_2013_2022...
std_err_mid_shelf_elon_2000_2018 ...
std_err_mid_shelf_vv_2000_2018 ...
std_err_shear_boundary_eshear_2000_2018 ...
std_err_pinning_point_elon_2013_2022...
std_err_pinning_point_vv_2000_2018 ...
std_err_pinning_point_elon_2000_2018 ...





%% -------------------------------------------------------------------------
% Cmobine the 2 datasets to get the total changes in the shear strain rates
% --------------------------------------------------------------------------
shear_boundary_eshear_time_series_2000_2023 = abs([(shear_boundary_eshear_time_series_2000_2018(1:12))...
    ;shear_boundary_eshear_time_series_2013_2022(2:end-1)]);

std_err_shear_boundary_eshear_2000_2023 = abs([(std_err_shear_boundary_eshear_2000_2018(1:12))...
    ;std_err_shear_boundary_eshear_2013_2022(2:end-1)]);

errorbar(shear_boundary_eshear_time_series_2000_2023,std_err_shear_boundary_eshear_2000_2023)
% 
% 
%% -----------------------
% mid-shelf time series
% -----------------------
mid_shelf_vv_time_series_2000_2023 = abs([mid_shelf_vv_time_series_2000_2018(2:13)...
    ;mid_shelf_vv_time_series_2013_2022(2:end-1)]);
% 
mid_shelf_elon_time_series_2000_2023 = ([mid_shelf_elon_time_series_2000_2018(2:13)...
    ;mid_shelf_elon_time_series_2013_2022(2:end-1)]);
% 
% 
% -----------------------------------------------
% standard error for the mid-shelf time series
% -----------------------------------------------
std_err_mid_shelf_vv_2000_2023 = abs([(std_err_mid_shelf_vv_2000_2018(2:13)) ...
    ;std_err_mid_shelf_vv_2013_2022(2:end-1)]);

std_err_mid_shelf_elon_2000_2023 = ([std_err_mid_shelf_elon_2000_2018(2:13)...
    ;(std_err_mid_shelf_elon_2013_2022(2:end-1))]);
% 
% 
% 
figure
bar(mid_shelf_elon_time_series_2000_2023)
hold on
errorbar(mid_shelf_elon_time_series_2000_2023,std_err_mid_shelf_elon_2000_2023)
% 
% figure
% errorbar(mid_shelf_vv_time_series_2000_2023,std_err_mid_shelf_vv_2000_2023)
% 
% 
% 
%% ---------------------------
% pinning point time series
% % ---------------------------
pinning_point_vv_time_series_2000_2023 = ([pinning_point_vv_time_series_2000_2018(2:13)...
    ;pinning_point_vv_time_series_2013_2022(2:end-1)]);

pinning_point_elon_time_series_2000_2023 = ([pinning_point_elon_time_series_2000_2018(2:13)...
    ;(pinning_point_elon_time_series_2013_2022(2:end-1))]);
% % 
% -------------------------------------------------
% standard error for the pinning point time series
% -------------------------------------------------
std_err_pinning_point_vv_2000_2023 = ([(std_err_pinning_point_vv_2000_2018(2:13)) ...
    ;std_err_pinning_point_vv_2013_2022(2:end-1)]);

std_err_pinning_point_elon_2000_2023 = ([std_err_pinning_point_elon_2000_2018(2:13)...
    ;(std_err_pinning_point_elon_2013_2022(2:end-1))]);


%
% figure
% errorbar(pinning_point_vv_time_series_2000_2023,std_err_pinning_point_vv_2000_2023)

figure
bar(pinning_point_elon_time_series_2000_2023)
hold on
errorbar(pinning_point_elon_time_series_2000_2023,std_err_pinning_point_elon_2000_2023)
% 
% close all

