%% User input
ispc= 0; %0 if linux, 1 if winOS
% linux work:
pathBase = '/run/media/reco/TOMO_GERD_2/APS/Data/Exp';
pathRest = '_D07T10N12_V0.068mlmin-1_C02.5E-4w%_multipleScanTest_dimax_2x_100mm_0.4msecExpTime_360DegPerSec_Rolling_100umLuAG_1mmC2mmGlass_2.657mrad_USArm1.1502_monoY_-10.0_AHutch/';
% mac home:
% stackBasePath = '/Volumes/Untitled/APS_Tomography/Exp041/Thres=0.0005,MedRad=0/bin_part_proj_0121_scan_'
% win OS:
% pathBase = 'f:\aps_2-bm_2016_08\exp';
% pathRest = '_d07t05n03_v0.068mlmin-1_c02.5e-4w%_multiplescantest_dimax_2x_100mm_0.4msecexptime_360degpersec_rolling_100umluag_1mmc2mmglass_2.657mrad_usarm1.1502_monoy_-10.0_ahutch\';
% ATTENTION: When changing experimentID also change the beginning of pathRest (sample id, V, c0)!
experimentIDt0 = 49;
experimentID = 50;
stackScanRange = 1:20; % actual number of time steps

% Which images from a single CT-scan (1...569) shall be read? 
imgRangeMin = 1;
imgRangeMax = 400;

% How many time steps may be skipped for a single trajectory
numSkip = 3

% How far can a particle travel in one time step?
maxLinkingDist = 50; %px

% particle volume, comes from particle size distribution plots and known
% particle diameter (d50 = 18µm = 3.6px)
minParticleSize = 7; %[px^3] equals d = 2.5px
maxParticleSize = 100; %[px^3] equals d = 5.8px


% to enable the drawing of found particle paths in an extra 3D image set
% following  variable to true.
boolDrawPathInImage = true;

%% Calculation
% iteration through TIME STEPS, for each time step, find all centers and
% establish one giant array containing all centers!


% initialize variable to store particle centers in
centerArray = cell(1,length(stackScanRange)) % as long as the number of frames we have!
k=1

for i=stackScanRange
    disp(['Time step is ' num2str(find(i==stackScanRange))]);
    disp('---');
    
    % read stack
    stackBaseFile = pathCreation(pathBase, pathRest, experimentID, i, ispc, 'part');
    stack = readtimeseries(stackBaseFile,'',[imgRangeMin imgRangeMax],0,0);
    % Read all "particles" from 3d image stack. 
    % Set minimum and maximum particle size to reduce computational cost
    % and minimize error.
    labels = label(stack>0,Inf,minParticleSize,maxParticleSize);

    % Next step: find particle position.
    msr = measure(labels,[],{'Center'},[],inf,0,0);
    
    try
        centerArray{k} = measures.center';
    catch
        centerArray{k} = [];
    end
    
    k = k+1;
    disp(['Number of particles found: ' num2str(size(measures.center))]);
    disp('=====================================================');
   
end


disp('Establishing link through all particle centers using Hungarian Linker');

[tracks,adjaceny_tracks] = simpletracker(centerArray,'MaxLinkingDist',maxLinkingDist,'MaxGapClosing',numSkip);

%rearrange the tracks so that it matches the style of particleTrackingGP:
% particles(i,j,k). i = 1:3 is x,y,z. j = 1:length(tracks) is the specific
% particle track and k = 1:length(stackScanRange) is the time step!
particles = zeros(3,length(tracks),length(stackScanRange));
for i = 1:length(tracks)
    for ii = 1:size(tracks{i},1)
        if(isnan(tracks{i}(ii))
            particles(:,i,ii) = 0
        else
            particles(:,i,ii) = centerArray{ii}(tracks{i}(ii),:);
        end
    end
end

   
    