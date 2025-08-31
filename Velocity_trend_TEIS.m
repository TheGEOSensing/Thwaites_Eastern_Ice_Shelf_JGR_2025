% Go to Karen's folder
cd('/Users/student/Desktop/ALL FILES/Datasets/Karen''s data on Thwaites Eastern Ice Shelf/')

info = geotiffinfo ('vv_500m_2000355_2002354.tif');

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


% %% Calculate the strain rates from the velocity files
% % ------------------------------------------------------- %
% for i = 1:size(vx,3)
%     [elon,etrans,eshear,eEff] = ...
%         strainRates( vx(:,:,i), vy(:,:,i), 2500);
%     % ------------------------------------------------------- %
%     % capture strain rates in individual 3d matrices
%     % ------------------------------------------------------- %
%     eps_long(:,:,i) = elon;
%     eps_trans(:,:,i) = etrans;
%     eps_shear(:,:,i) = eshear;
%     eps_Eff(:,:,i) = eEff;    
% end
% 
% 
% flow_angle = atan2d(vx, vy);
% 
% clearvars -except eps_long eps_trans eps_Eff eps_shear...
%     vx vy vv R %(for strain rates)



%% Let's see the trend in velocity in different times

% First, set up the different periods according to the literature
% ------------------------------------------------------- %
TWIT_Acceleration_velocity = vv(:,:,2:6); % 2002 - 2006
Shear_break_velocity = vv(:,:,7:11);      % 2007 - 2011
Post_shear_break_velocity = vv(:,:,11:16);   % 2011- 2016


% now calculate the trend from the velocity data in different times
% ------------------------------------------------------- %
TWIT_Acc_velocity_trend = calculateTrend3D (TWIT_Acceleration_velocity);
Shear_break_velocity_trend = calculateTrend3D (Shear_break_velocity);
Post_shear_break_velocity_trend = calculateTrend3D (Post_shear_break_velocity);



% Plot the trend in veloity one by one
% ------------------------------------------------------- %
figure(1)
subplot(1,3,1)
imagesc(flipud(TWIT_Acc_velocity_trend)*133.4);
colormap(brewermap(20, 'RdBu'));
caxis([-50 50])


subplot(1,3,2)
imagesc(flipud(Shear_break_velocity_trend)*133.4);
colormap(brewermap(20, 'RdBu'));
caxis([-50 50])


subplot(1,3,3)
imagesc(flipud(Post_shear_break_velocity_trend)*133.4);
colormap(brewermap(20, 'RdBu')); 
caxis([-50 50])
colorbar

%% export the geotiff files of trend of velocity
% ------------------------------------------------------- %
cd('/Users/student/Desktop/Objective 1/Trend of velocity geotiffs/')

geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
tiffTags = struct('TileLength',1024,'TileWidth',1024);

outfilename = 'Post_shear_break_velocity_trend.tif';
geotiffwrite(outfilename,flipud(Post_shear_break_velocity_trend *133.4),R,'TiffType','bigtiff', ...
                             'GeoKeyDirectoryTag',geoTags, ...
                             'TiffTags',tiffTags)


%% Trend for the last phase
% --------------------------
cd('/Users/student/Desktop/Objective 1/Its_live vleocity_image_pairs')
load('Thwaites_velocity_image_Pair_products.mat')

velocity_2017_2022 = velocity(:,:,4:10);
velocity_trend_2017_2022 = calculateTrend3D(velocity_2017_2022);



cd '/Users/student/Desktop/Objective 1/Trend of velocity geotiffs'/
geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
tiffTags = struct('TileLength',1024,'TileWidth',1024);

outfilename = ['Velocity_Trend_2017_2022.tif'];
geotiffwrite(outfilename, velocity_trend_2017_2022, R_vel, 'TiffType', 'bigtiff', ...
                 'GeoKeyDirectoryTag', geoTags, 'TiffTags', tiffTags);


%%

velocity_2013_2017 = velocity(:,:,1:4);
velocity_trend_2013_2017 = calculateTrend3D(velocity_2013_2017);

figure
imagesc((velocity_trend_2013_2017));
colormap(brewermap(20, 'RdBu')); 
caxis([-50 50])
colorbar


cd '/Users/student/Desktop/Objective 1/Trend of velocity geotiffs'/
geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
tiffTags = struct('TileLength',1024,'TileWidth',1024);

outfilename = ['Velocity_Trend_2013_2017.tif'];
geotiffwrite(outfilename, velocity_trend_2013_2017, R_vel, 'TiffType', 'bigtiff', ...
                 'GeoKeyDirectoryTag', geoTags, 'TiffTags', tiffTags);

