
% Script opens a series of stacks for particle analysis
% Aim is to test the binarisation parameters for particle recognition
% (aps_2-bm_2016_08 experiments). 
% Uses DIPImage library to recognize particles and return number of
% recognized particles and a vector of particle sizes.
% number of particles is a vector of size (X,1) with X being the number
% of scanned stackes. sizeDist is particle size distribution (Y,X),
% Y is the number of particle size classes (length of vector
% sizeDistClasses)

% User input
%linux work:
stackBasePath = '/home/gp/Pesch/Ph.D./unshared/APS_Tomographie/Exp041/Thres=0.0005,MedRad=2/bin_part_proj_0121_scan_';
%mac home:
%stackBasePath = '/Volumes/Untitled/APS_Tomography/Exp041/Thres=0.0005,MedRad=0/bin_part_proj_0121_scan_'
project = '0121'
stackScanRange = [1:5]


sizeMin=1;
sizeMax=50;
numSizes=50;

% Some preliminiary calculation
sizeDistClasses = linspace(sizeMin,sizeMax,numSizes);
sizeDist = zeros(numSizes,length(stackScanRange));


z = 1; %index for particle size classes!
for i=stackScanRange
    % path creation is a little holprig because I made it quick and dirty:
    stackBaseFile = [stackBasePath, num2str(i,'%02.0f'), '/bin_part_proj_', project, '_scan_', num2str(i,'%02.0f'), '_1*.tif'];
    %stack = readtimeseries(stackBaseFile,'',[],0,1);
    stack = readtimeseries(stackBaseFile,'',[0 183],0,1);

    % read all "particles" from 3d image stack. no restrictions on
    % dimension or anything else.
    % no false positive removal at this stage!
    labels = label(stack>0,Inf,9,11);

    % next step is to apply some measurements to the identified
    % "particles":
    msr = measure(labels,[],{'Size'},[],inf,0,0);

    numberParticles=max(msr.id);
    for ii=1:numberParticles
        sizeClassIndex = find(sizeDistClasses>msr.size(ii),1,'first');
        if (~isempty(sizeClassIndex))
            sizeDist(sizeClassIndex,z) = sizeDist(sizeClassIndex,z)+1;
        end  
    end
    figure;
    stem(sizeDistClasses,sizeDist(:,z));
z=z+1;    
end

