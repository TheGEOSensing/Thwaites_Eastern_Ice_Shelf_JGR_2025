% Go to Karen's folder
cd('/Users/student/Desktop/ALL FILES/Datasets/Karen''s data on Thwaites Eastern Ice Shelf/')

% Read all the velocity images one by ne using a for loop
vv_files = dir('vv_*.tif');
vx_files = dir('vx_*.tif');
vy_files = dir('vy_*.tif');

for i = 1:size(vv_files)

% the vv files as velocity image
% -------------------------------- %
    image1 = vv_files(i).name;
    image = readgeoraster(image1);
    vv(:,:,i) = flipud(image);

% vx files as the longitudinal velocity
% ------------------------------------- %
    image1 = vx_files(i).name;
    image = readgeoraster(image1);
    vx(:,:,i) = flipud(image);

% vy files as the transverse velocity
% -------------------------------------- %
    image1 = vy_files(i).name;
    [image,R] = readgeoraster(image1);
    vy(:,:,i) = flipud(image);
end

info = geotiffinfo ('vv_500m_2000355_2002354.tif');


%% Calculate the strain rates from the velocity files
% ------------------------------------------------------- %
for i = 1:size(vx,3)
    [elon,etrans,eshear,exGrid,eyGrid,exyGrid,eEff] = ...
        strainRates( vx(:,:,i), vy(:,:,i), 500);
    % ------------------------------------------------------- %
    % capture strain rates in individual 3d matrices
    % ------------------------------------------------------- %
    eps_long(:,:,i) = elon;
    eps_trans(:,:,i) = etrans;
    eps_shear(:,:,i) = eshear;
    eps_Eff(:,:,i) = eEff;    

    eps_xx(:,:,i) = exGrid;
    eps_yy(:,:,i) = eyGrid;
    eps_xy(:,:,i) = exyGrid;    

end

divergence = (eps_xx + eps_yy);

clearvars -except eps_long eps_trans eps_Eff eps_shear...
    vx vy vv R divergence %(for strain rates)



%% export the geotiff files of trend of velocity
% ------------------------------------------------------- %
%cd('/Users/student/Desktop/Objective 1/Shear_strain_rate_Geotiffs')
%cd '/Users/student/Desktop/Objective 1/TEIS_geotiffs_2002_2020_Karen''s_data/Shear_Strain_Rates'
% 
% geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
% tiffTags = struct('TileLength',1024,'TileWidth',1024);
% 
% for i = 1:size(eps_shear, 3)
%     outfilename = ['Divergence_' num2str(2001 + i) '.tif'];
%     geotiffwrite(outfilename, flipud(divergence(:, :, i)), R, 'TiffType', 'bigtiff', ...
%                  'GeoKeyDirectoryTag', geoTags, 'TiffTags', tiffTags);
% end


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

clearvars -except eps_long eps_trans eps_Eff eps_shear...
    vx vy vv R X Y Latitude Longitude  divergence %(for strain rates)


%% GET THE SHEAR BOUNDARY SHAPEFILE
% ------------------------------------------------------- %

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


for i = 1:size(vv,3)
    % --------------------------------------------------- %
    % make a masked for each of the datasets
    % multiply the masked with the matrix to clip
    % --------------------------------------------------- %
    masked_vv = (vv(:,:,i)).* mask;        masked_vv(masked_vv == 0) = NaN;
    masked_vx = (vx(:,:,i)).* mask;        masked_vx(masked_vx == 0) = NaN;
    masked_vy = (vy(:,:,i)).* mask;        masked_vy(masked_vy == 0) = NaN;

    masked_elon = (eps_long(:,:,i)).* mask;       
    masked_ehsear = (eps_shear(:,:,i)).* mask;    
    masked_etrans = (eps_trans(:,:,i)).* mask;   

    % --------------------------------------------------- %
    % combine the temporal datasets together
    % --------------------------------------------------- %
    shear_boundary_vv(:,:,i) = masked_vv;
    shear_boundary_vx(:,:,i) = masked_vx;
    shear_boundary_vy(:,:,i) = masked_vy;

    shear_boundary_elon(:,:,i) = masked_elon;
    shear_boundary_eshear(:,:,i) = masked_ehsear;
    shear_boundary_etrans(:,:,i) = masked_etrans;
end

figure
pcolorps(Latitude,Longitude, (masked_ehsear(:,:,1)));

figure
pcolorps(Latitude,Longitude, (abs(eps_shear(:,:,1))));

%% LET'S DO IT FOR THE MID-SHELF AREA (LONGITUDINAL STRAIN RATE)
% --------------------------------------------------------------- %

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

for i = 1:size(vv,3)
    % --------------------------------------------------- %
    % make a masked for each of the datasets
    % multiply the masked with the matrix to clip
    % --------------------------------------------------- %
    masked_vv = (vv(:,:,i)).* mask;        masked_vv(masked_vv == 0) = NaN;
    masked_vx = (vx(:,:,i)).* mask;        masked_vx(masked_vx == 0) = NaN;
    masked_vy = (vy(:,:,i)).* mask;        masked_vy(masked_vy == 0) = NaN;

    masked_elon = (eps_long(:,:,i)).* mask;    
    masked_ehsear = (eps_shear(:,:,i)).* mask;   
    masked_etrans = (eps_trans(:,:,i)).* mask;   

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

figure
pcolorps(Latitude,Longitude, (masked_ehsear(:,:,1)));

figure
pcolorps(Latitude,Longitude, (abs(eps_shear(:,:,1))));


%% LET'S DO IT FOR THE MID-SHELF AREA (LONGITUDINAL STRAIN RATE)
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

for i = 1:size(vv,3)
    % --------------------------------------------------- %
    % make a masked for each of the datasets
    % multiply the masked with the matrix to clip
    % --------------------------------------------------- %
    masked_vv = (vv(:,:,i)).* mask;        masked_vv(masked_vv == 0) = NaN;
    masked_vx = (vx(:,:,i)).* mask;        masked_vx(masked_vx == 0) = NaN;
    masked_vy = (vy(:,:,i)).* mask;        masked_vy(masked_vy == 0) = NaN;

    masked_elon = (eps_long(:,:,i)).* mask;    
    masked_ehsear = (eps_shear(:,:,i)).* mask;   
    masked_etrans = (eps_trans(:,:,i)).* mask;   

    % --------------------------------------------------- %
    % combine the temporal datasets together
    % --------------------------------------------------- %
    pinning_point_vv(:,:,i) = masked_vv;

    pinning_point_elon(:,:,i) = masked_elon;
    pinning_point_eshear(:,:,i) = masked_ehsear;
    pinning_point_etrans(:,:,i) = masked_etrans;
end


figure
pcolorps(Latitude,Longitude, (masked_ehsear(:,:,1)));


figure
pcolorps(Latitude,Longitude, (abs(eps_shear(:,:,1))));

%% ------------------------------------------------------- %
% Let's plot a figure
% ------------------------------------------------------- %

subplot(1,2,1)
pcolorps(Latitude, Longitude, shear_boundary_eshear(:,:,12))
colormap((brewermap([],"RdBu")));
caxis([-0.2*10^-4  0.2*10^-4])
%caxis([0 15])
%plotps(ry, rx, 'k', 'linewidth', 3.5);
antbounds('gl','polyshape','facealpha',1, 'facecolor',[.8 .8 .8],'linewidth',1.5)
antbounds('coast','black','linewidth',1.5)
set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

%
subplot(1,2,2)
pcolorps(Latitude, Longitude, eps_shear(:,:,12))
colormap((brewermap([],"RdBu")));
caxis([-2*10^-4  2*10^-4])
%caxis([0 15])
%plotps(ry, rx, 'k', 'linewidth', 3.5);
antbounds('gl','polyshape','facealpha',1, 'facecolor',[.8 .8 .8],'linewidth',1.5)
antbounds('coast','black','linewidth',1.5)
set(gca, 'xtick', []);
set(gca, 'ytick', []);
box on;

%% ------------------------------------------------------- %
% Clear the workspace as much as possible
% ------------------------------------------------------- %

clearvars -except ...
eps_long eps_trans eps_Eff eps_shear...
vx vy vv R ...
shear_boundary_eshear shear_boundary_elon shear_boundary_etrans ...
shear_boundary_vv shear_boundary_vx shear_boundary_vy ...
Latitude Longitude X Y rx ry ...
mid_shelf_etrans  mid_shelf_eshear  mid_shelf_elon ... 
mid_shelf_vv  mid_shelf_vx  mid_shelf_vy ...
pinning_point_vv pinning_point_elon pinning_point_etrans pinning_point_eshear


%% Let's try the changes in shear strain rates 
% ------------------------------------------------------- %

% ------------------------------------------------------- %
% shear zone time series of strain rates
% ------------------------------------------------------- %
shear_boundary_eshear_time_series_2000_2018 = squeeze(nansum(nansum(abs(shear_boundary_eshear),1),2));

% ------------------------------------------------------- %
% Mid shelf time series of velocity
% ------------------------------------------------------- %
mid_shelf_vv_time_series_2000_2018 = squeeze(nanmean(nanmean(mid_shelf_vv,1),2));


% ------------------------------------------------------- %
% Mid shelf time series of strain rates
% ------------------------------------------------------- %
mid_shelf_elon_time_series_2000_2018 = squeeze(nansum(nansum(mid_shelf_elon,1),2));


% ------------------------------------------------------- %
% Pinning Point time series of velocity
% ------------------------------------------------------- %
pinning_point_vv_time_series_2000_2018 = squeeze(nanmean(nanmean(pinning_point_vv,1),2));


% ------------------------------------------------------- %
% Pinning Point time series of strain rates
% ------------------------------------------------------- %
pinning_point_elon_time_series_2000_2018 = squeeze(nansum(nansum(pinning_point_elon,1),2));



%% ------------------------------------------------------- %
% Standard Deviation for the calculation of the errorbar 
% -------------------------------------------------------- %
% number of sample calculation %
% -------------------------------
n_shear_zone = sum(~isnan(shear_boundary_vv(:)))/size(shear_boundary_vv,3);
n_mid_shelf = sum(~isnan(mid_shelf_vv(:)))/size(shear_boundary_vv,3);
n_pinning_point = sum(~isnan(pinning_point_vv(:)))/size(shear_boundary_vv,3);





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
% ---------------
std_err_mid_shelf_elon_2000_2018 = squeeze(nanstd(mid_shelf_elon, [], 1:2)/root_n_1_mid_shelf) * n_mid_shelf;
std_err_mid_shelf_vv_2000_2018 = squeeze(nanstd(mid_shelf_vv, [], 1:2)/root_n_1_mid_shelf);

% shear boundary
% ---------------
std_err_shear_boundary_eshear_2000_2018 = squeeze(nanstd(abs(shear_boundary_eshear), [], 1:2)/root_n_1_shear_zone) * n_shear_zone;
 
% pinning point
% ---------------
std_err_pinning_point_vv_2000_2018 = squeeze(nanstd(pinning_point_vv, [], 1:2)/root_n_1_pinning_point);
std_err_pinning_point_elon_2000_2018 = squeeze(nanstd(pinning_point_elon, [], 1:2)/root_n_1_pinning_point) * n_pinning_point;



%%
clearvars -except ...
mid_shelf_etrans_time_series_2000_2018...
mid_shelf_elon_time_series_2000_2018...
mid_shelf_eshear_time_series_2000_2018...
mid_shelf_vv_time_series_2000_2018 ...
shear_boundary_etrans_time_series_2000_2018...
shear_boundary_elon_time_series_2000_2018...
shear_boundary_eshear_time_series_2000_2018...
pinning_point_eshear_time_series_2000_2018 ...
pinning_point_etrans_time_series_2000_2018 ...
pinning_point_elon_time_series_2000_2018 ...
pinning_point_vv_time_series_2000_2018 ...
std_err_mid_shelf_elon_2000_2018 ...
std_err_mid_shelf_vv_2000_2018 ...
std_err_shear_boundary_eshear_2000_2018 ...
std_err_pinning_point_vv_2000_2018 ...
std_err_pinning_point_elon_2000_2018 ...

%%
close all

%%
figure
errorbar(shear_boundary_eshear_time_series_2000_2018,std_err_shear_boundary_eshear_2000_2018)
title('Shear Zone eShear')

figure
errorbar(mid_shelf_elon_time_series_2000_2018,std_err_mid_shelf_elon_2000_2018)
title(' Mid Shelf eLong')

figure
errorbar(pinning_point_elon_time_series_2000_2018,std_err_pinning_point_elon_2000_2018)
title(' Pinning Point eLong')

