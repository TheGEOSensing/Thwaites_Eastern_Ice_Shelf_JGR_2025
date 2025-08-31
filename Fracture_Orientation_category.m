
% Location of the fracture shapefiles
cd '/Users/rishi/Desktop/Objective 1/Clipped_fractures_TEIS/Fracture_Orientation_epsg3031'

files = {'fractures_2002.shp',...
    'fractures_2006.shp', ...
    'fractures_2010.shp',...
    'fractures_2014.shp',...
    'fractures_2018.shp',...
    'fractures_2021.shp'};


% Initialize cell arrays to store Orientation values for each shapefile
Orientation_grt_10 = cell(1, numel(files));
Orientation_8_to_10 = cell(1, numel(files));
Orientation_6_to_8 = cell(1, numel(files));
Orientation_4_to_6 = cell(1, numel(files));
Orientation_2_to_4 = cell(1, numel(files));
Orientation_less_2 = cell(1, numel(files));


% Loop through each shapefile
for fileIdx = 1:numel(files)
    % Read the shapefile and get the fracture orientation from it
    fracture = shaperead(files{fileIdx});
    
    % Extract Orientation based on different length ranges for the current shapefile
    for i = 1:numel(fracture)
        if fracture(i).length > 10
            Orientation_grt_10{fileIdx} = [Orientation_grt_10{fileIdx}; fracture(i).Orientatio];
        elseif fracture(i).length >= 8 && fracture(i).length < 10
            Orientation_8_to_10{fileIdx} = [Orientation_8_to_10{fileIdx}; fracture(i).Orientatio];
        elseif fracture(i).length >= 6 && fracture(i).length < 8
            Orientation_6_to_8{fileIdx} = [Orientation_6_to_8{fileIdx}; fracture(i).Orientatio];
        elseif fracture(i).length >= 4 && fracture(i).length < 6
            Orientation_4_to_6{fileIdx} = [Orientation_4_to_6{fileIdx}; fracture(i).Orientatio]; 
        elseif fracture(i).length >= 2 && fracture(i).length < 4
            Orientation_2_to_4{fileIdx} = [Orientation_2_to_4{fileIdx}; fracture(i).Orientatio];
        elseif fracture(i).length < 2
            Orientation_less_2{fileIdx} = [Orientation_less_2{fileIdx}; fracture(i).Orientatio];        
        end
    end
end


% Seperate the years based on the cells created in the last code 

% Define the suffix for variable names
years = {'2002', '2006', '2010', '2014', '2018', '2021'};

for fileIdx = 1:numel(files)
    
    varName = ['Orientation_' years{fileIdx} '_grt_10']; % Define variable name
    eval([varName ' = vertcat(Orientation_grt_10{' num2str(fileIdx) '});']); % Concatenate cell arrays

    varName = ['Orientation_' years{fileIdx} '_8_10'];
    eval([varName ' = vertcat(Orientation_8_to_10{' num2str(fileIdx) '});']);

    varName = ['Orientation_' years{fileIdx} '_6_8'];
    eval([varName ' = vertcat(Orientation_6_to_8{' num2str(fileIdx) '});']);

    varName = ['Orientation_' years{fileIdx} '_4_6'];
    eval([varName ' = vertcat(Orientation_4_to_6{' num2str(fileIdx) '});']);

    varName = ['Orientation_' years{fileIdx} '_2_4'];
    eval([varName ' = vertcat(Orientation_2_to_4{' num2str(fileIdx) '});']);

    varName = ['Orientation_' years{fileIdx} '_less_2'];
    eval([varName ' = vertcat(Orientation_less_2{' num2str(fileIdx) '});']);

end


%% polar histogram plot for the fracture orientation of fracture length <4km
% ----------------------------------------------------------------------------

Orientation_2002_less_4 = [Orientation_2002_less_2; Orientation_2002_2_4];
Orientation_2006_less_4 = [Orientation_2006_less_2; Orientation_2006_2_4];
Orientation_2010_less_4 = [Orientation_2010_less_2; Orientation_2010_2_4];
Orientation_2014_less_4 = [Orientation_2014_less_2; Orientation_2014_2_4];
Orientation_2018_less_4 = [Orientation_2018_less_2; Orientation_2018_2_4];
Orientation_2021_less_4 = [Orientation_2021_less_2; Orientation_2021_2_4];

% List of orientations for different years
Orientations = {Orientation_2002_less_4, Orientation_2006_less_4, Orientation_2010_less_4, ...
    Orientation_2014_less_4, Orientation_2018_less_4 , Orientation_2021_less_4};


Orientation_years = [2002, 2006, 2010, 2014, 2018, 2021];

% Create a figure with subplots for each orientation
% --------------------------------------------------
f = figure;
f.Position = [100 100 1200 800];

for i = 1:numel(Orientations)

    % design the polar axis for the graphics
    % ---------------------------------------
    axesHandle(i) = subplot(2, 3, i);
    polarAxesHandle(i) = polaraxes('Units', axesHandle(i).Units, 'Position', axesHandle(i).Position);
    delete(axesHandle(i));
    
    % plot the orientation values on the polar histogram axis
    % -------------------------------------------------------
    p3 = polarhistogram(polarAxesHandle(i), deg2rad(mod(Orientations{i},180)), 12,...
        'FaceAlpha', 0.5,'BinLimits', [0, pi]);
    p3.LineWidth = 1; 
    polarAxesHandle(i).FontSize = 12;
    polarAxesHandle(i).FontWeight = 'bold';
    polarAxesHandle(i).GridAlpha = 0.4;
    polarAxesHandle(i).ThetaDir = 'clockwise';
    polarAxesHandle(i).ThetaZeroLocation = 'top';
    rlim([0 30]);
    rticks([ 10 15 20 25 30]);

    % Add title based on the year
    % ----------------------------
    title(polarAxesHandle(i), num2str(Orientation_years(i)), 'FontWeight', 'bold', 'FontSize',18);
    subtitle(polarAxesHandle(i), 'Fracture Orientation (\leq4km)', 'FontWeight', 'bold', 'FontSize',14);
end
