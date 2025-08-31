% Smooth GPS velocities with 10-day median filter
speedSmooth = smoothdata(cavccont.Speed,'movmedian',961);

% Read ITS_LIVE average geotiff
[vels,vel_info] = readgeoraster('/Volumes/Linus/ItsLive/Antarctic_mosaics/ANT_G0240_0000_vv.tif');
[vvcrop,vvcrop_info] = mapcrop(vels,vel_info,[min(cavccont.PolarStereoX)-1200,max(cavccont.PolarStereoX)+1200],[min(cavccont.PolarStereoY)-1200,max(cavccont.PolarStereoY)+1200]);

% Plot points over map background
figure(1)
[xcrop,ycrop] = worldGrid(vvcrop_info);
h = pcolor(xcrop,ycrop,vvcrop);
hold on
set(h, 'EdgeColor','none');
scatter(cavccont.PolarStereoX,cavccont.PolarStereoY)


% Extract background velocities at GPS points
avVels = mapinterp(vvcrop,vvcrop_info,cavccont.PolarStereoX,cavccont.PolarStereoY);
avVels = avVels/365.25;
diffAvVels = diff(avVels);
VelCorr = ones(size(avVels));
VelCorr(1) = speedSmooth(1);

% Calculate cumulative impact of spatial gradient along path for each GPS
% point and subtract it from original speeds
corrValue = 0;
for i = 2:length(speedSmooth)
    corrValue = corrValue + diffAvVels(i-1);
    VelCorr(i) = speedSmooth(i) - corrValue;
end

% Plot uncorrected and corrected speeds over time
figure(2)
plot(cavccont.DOY2020,speedSmooth)
hold on
plot(cavccont.DOY2020,VelCorr)

