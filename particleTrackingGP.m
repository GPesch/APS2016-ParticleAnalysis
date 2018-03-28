%% changelog.
% original version  particleTrackingGP
% version 0_0_1     added subfunction calcLengthAndTimeSteps() for further evaluation
%                   - truncate 'particles' array to particles which were
%                   found in e.g. 10 time steps (-> user input)
%                   - calculate the length of the path as additional info
%                   - meta data is automatically saved in folder
% version 0_0_2     added subfunction drawPathInImage() (-> user input)
%                   -> needs calcLengthAndTimeSteps() becasue of idsOfInterest

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

% specify the search radius 
r0 = 15; %[px]
dr = 5; %[px]
rmax = 50; %[px] equals (rmax-r0)/dr = 7 + 1 iterations


% particle volume, comes from particle size distribution plots and known
% particle diameter (d50 = 18µm = 3.6px)
minParticleSize = 7; %[px^3] equals d = 2.5px
maxParticleSize = 100; %[px^3] equals d = 5.8px

% relevant time steps. The user is interested in particles that are
% found/matched in that many successive time steps.
% ATTENTION: must be <= size(stackScanRange,2) (i.e. number of time steps)
% Do this calculation (true/false)?
boolRelevTimeSteps = true;
relevTimeSteps = 10;

% to enable the drawing of found particle paths in an extra 3D image set
% following  variable to true.
boolDrawPathInImage = true;

%% Calculation
% iteration through TIME STEPS
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
    % Sort measure array (if there are more than 1008 IDs, msr.id(i)
    % [=entry i in array] is not equal to msr(i).id [=ID==i] - also true
    % for center, etc.).
    [msrIDs, msrCenters]=treatArray(msr,'sort',0);
    
    numberParticles=max(msr.id);
    
    % At the first image we fill an array with initial particle positions
    % we then try to find in the other time steps to find partners.
    
    % first time step t=t1?
    if i == stackScanRange(1)
        % initialize particle array
        % three-dimensional; dim 1 = [x y z], dim 2 = particle number ii,
        % dim 3 = time step i
        %
        % t1:                       t2:
        %    |p1 |p2 |p3 |p4 ...       |p1 |p2 |p3 |p4 ...
        %    ----------------          ----------------
        %  x |x1 |x2 |x3 |x4 ...     x |x1 |x2 |x3 |x4 ...
        %  y |y1 |y2 |y3 |y4 ...     y |y1 |y2 |y3 |y4 ...
        %  z |z1 |z2 |z3 |z4 ...     z |z1 |z2 |z3 |z4 ...
        %
        particles = zeros(3,numberParticles,max(stackScanRange));
        particles(:,:,1) = msrCenters(:,:);
    else % t>t_1
        % start off with small radius, find _SINGLE_ partners. 
        % populate three arrays: found, skip and multipleFound. 
        %   found           contains all particles from time step i-1 
        %                   that found exactly one partner in time step i
        %                   and vice versa
        %                   o --> o     (good)
        %   skip            contains all particles from time step i-1 that 
        %                   found multiple partners in time step i
        %                   o ---> o o  (bad)
        %   multipleFound   contains all particles from time step i-1 where
        %                   at least two particles from time step i-1 
        %                   found the same particle in time step i
        %                   o o ---> o  (bad)
        % Particles that do not appear in one of those three arrays had no
        % match in time step i. The radius is increased to deal with those.
        
        found = int16.empty;
        skip = int16.empty;
        multipleFound = int16.empty;
            
        % Particles Buffer is an auxiliary array that stores the particle
        % IDs of the current time step in order to check whether there was a single match.
        % This is true for case 'found' and 'multipleFound'. A
        % differentiation between both cases is done after all particles
        % have been checked within one radius-loop.
        %                      
        %    |p1 |p2 |p3 |p4 ...
        %    ----------------  
        %  id|id1|id2|id3|id4...
        %  x |x1 |x2 |x3 |x4 ...
        %  y |y1 |y2 |y3 |y4 ...
        %  z |z1 |z2 |z3 |z4 ...
        %
        particlesBuffer=zeros(4,size(particles,2));
        
        for ri=r0:dr:rmax
            disp(['Search radius ri =               : ' num2str(ri) ]);
            disp(['Considering                      : ' num2str(length(nonzeros(particles(1,:,i-1)))) ' / ' num2str(size(particles,2)) ' particles' ])
            disp(['Possible partners                : ' num2str(length(msr.Id)) ]);
            
            for ii=1:size(particles,2) % iterate through all particles found in t=t1
                % skip all particles that already have partners (1 or
                % more). E.g. in a previous time step particle ii
                % was found to have several possible matches. Hence its id
                % was added to 'skip'. In the current time step t this
                % particle is not further considered due to the following if-statement.
                if ismember(ii,horzcat(found,skip,multipleFound))
                    continue
                end
                
                % skip all particles that didn't have a partner the
                % timestep before.
                if sum(particles(:,ii,i-1) == 0) %% any before
                    continue
                end
                
                % EXAMPLE:
                %
                % The position of particle ii=5 in time step t=t1 (i=1) is 
                % known and stored in particles(:,5,1). In time step t=t2 
                % (i=2) a single particle is found in the vicinity of that
                % particle. Its position is stored in particles(:,5,2). In
                % time step t=t3, the position of t2 is the basis to find
                % another match and so on.
                
                [nearbyIds,nearbyCenters] = findNearbyParticles(particles(:,ii,i-1),msrCenters,ri);
                % particles(:,ii,i-1) are all three dimensions of particle ii at the last time step (i-1).
                
                % if length(nearbyIds)==0, there was no nearby particle. No
                % action is taken in this case. If ri is increased, a
                % nearby particle might be found.
                
                if length(nearbyIds)==1
                    % algorithm found precisely one particle
                    % The position of particlesBuffer equals old particle id (t_{i-1}), the
                    % value equals nearbyIds (=new particle id, t_i).
                    particlesBuffer(1,ii)=nearbyIds;
                    particlesBuffer(2:4,ii)=nearbyCenters;
                elseif length(nearbyIds) > 1
                    % algorithm found more than one particle partner,
                    % skip this particle
                    skip = horzcat(skip,ii);
                end
            end
                            
            % After all particles (time step i-1) have been checked,
            % clarify which particles from time step i are matched once/multiple times.
            for kk=1:max(particlesBuffer(1,:)) % max nearbyIds in particlesBuffer
                % The position (row 1, col xy) of particlesBuffer equals ii (particle id of
                % the first time step; here tmp), the value equals nearbyIds 
                % (particle ids of the curretn time step; here kk).
                tmp = find(particlesBuffer(1,:)==kk);
                if size(tmp,2)==1 % one particle has been matched once
                    if ~ismember(tmp,found) % if it is not already in found
                        found = horzcat(found, tmp); % add it
                    end
                elseif size(tmp,2)>1 % one particle has been matched multiple times
                    for bb=1:size(tmp,2)
                        if ~ismember(tmp(bb),multipleFound) % if it is not already in multipleFound
                            multipleFound = horzcat(multipleFound,tmp(bb)); % add it
                        end
                    end
                end
            end
            % Sort found.
            if(~isempty(found))
                found = (sortrows(found',1))';
                % Remove it from measure array by putting zeros in the
                % array at the specific row. In the next iteration of the
                % radius those are not considered any more.
                [~,msrCenters]=treatArray(msrCenters,'remove',particlesBuffer(1,found));
            end
            
            disp(['Found partners in ri  o   -> o   : ' num2str(size(found,2)) ]);
            disp(['Too many partners     o   -> o o : ' num2str(size(skip,2)) ]);
            disp(['Found same particle   o o -> o   : ' num2str(size(multipleFound,2)) ]);
            disp(['No match found                   : ' num2str(length(nonzeros(particles(1,:,i-1)))-size(multipleFound,2)-size(skip,2)-size(found,2)) ]);
            disp('---');
        end
        % Save new particle position.
        particles(:,found,i) = particlesBuffer(2:4,found);
        disp('=====================================================');
    end     
end

%save 'particles' array
dlmwrite([pathBase num2str(experimentID,'%03.0f') pathRest 'particlesTracked_particleTrackingGP_0_0_1_FiltEval0_2_1.txt'],particles,'delimiter','\t','precision',8)

%% Check if individual particle displacement vectors penetrate filter wall
% Particle displacement vectors must not penetrate the filter wall. If so,
% the particle tracking algorithm found (most likely) an erroneous match.

% Read binarized image of porous filter (1==bulk, 0==pore).
stackBaseFileT0 = pathCreation(pathBase, pathRest, experimentIDt0, 1, ispc, 'filt');
stackFilt = readtimeseries(stackBaseFileT0,'',[imgRangeMin imgRangeMax],0,0)>0;

% Loop through particles
for bb=1:size(particles,2)
    x = reshape(particles(1,bb,:),[length(stackScanRange) 1]);
    y = reshape(particles(2,bb,:),[length(stackScanRange) 1]);
    z = reshape(particles(3,bb,:),[length(stackScanRange) 1]);
    
    if any(x + y + z == 0) % is not in 'found'
        continue
    end

    for cc=1:length(stackScanRange)-1
        % vector of particle bb between two time steps (cc and cc+1)
        xx = [x(cc+1)-x(cc), y(cc+1)-y(cc), z(cc+1)-z(cc)];
        % length (norm) of xx (rounded down)
        normVec = floor(norm(xx));
        if normVec<1    % particle displacement is smaller than a single px (5mu)...skip that
            continue
        end
        dirVec = xx/normVec;
        
        % Check if vector xx penetrates filter wall. To do so check
        % pixel-wise along the displacement vector if particle leaves a
        % pore (penetrates the filter).
        for dd=1:normVec
            point = round([x(cc), y(cc), z(cc)]+dd*dirVec);
            if dip_array(stackFilt(point(1),point(2),point(3)))==1 % == filter
                % If so, put flag in particles(:,bb,cc). x, y and z will be
                % negative (out of image bounds) which can be used later on in display section.
                particles(:,bb,cc) = -particles(:,bb,cc);
                break
            end
        end
    end
end

% save 'particles' array
% ATTENTION: saved array has only 2 dimensions. (3, particles x
% time_steps). Must be rearranged if reloaded (using particles=dlmread(path)).
dlmwrite([pathBase num2str(experimentID,'%03.0f') pathRest 'particlesTracked_particleTrackingGP_0_0_1_FiltEval0_2_1.txt'],particles,'delimiter','\t','precision',8);

%% Further evaluation of 'particles' array (e.g. length of path)

if boolRelevTimeSteps % -> user input
    % do calculations. Results are saved as txt-files. Two results are
    % returned (the particle ids that match the conditions - i.e. those 
    % were found in x time steps - and the corresponding path lengths).
    [idsOfInterest, pathLength] = calcLengthAndTimeSteps(particles, relevTimeSteps, stackScanRange, pathBase, pathRest, experimentID);
end
% look at path length. e.g.: hist(pathLength(:,3),[1 2 5 10 20 50 100 200 500])


if (boolDrawPathInImage && boolRelevTimeSteps) % -> user input
    % draw found particle paths in a 3D image. Main input value is
    % 'particle' array (without negative values -> abs(particles)).
    % No return value.
    drawPathInImage(abs(particles), stackScanRange, idsOfInterest, (imgRangeMax - imgRangeMin + 1), pathBase, pathRest, experimentID);
end

%% Display
% strip array of all particles that contain 0 coordinates:
figure;
hold on
j = 0;
for aa=1:size(particles,2)
    % rearrange
    x = reshape(particles(1,aa,:),[length(stackScanRange) 1]);
    y = reshape(particles(2,aa,:),[length(stackScanRange) 1]);
    z = reshape(particles(3,aa,:),[length(stackScanRange) 1]);

    if any(x + y + z == 0) % is not in 'found'
        continue
    elseif any (x + y + z < 0) % at some point penetrates Filter
        j = j+1;
        hitFilt = find(x + y + z < 0); % in which of the displacement-vector sections is the filter hit (penetrated)
        % Make values positive again for proper plotting.
        x(hitFilt) = -x(hitFilt);
        y(hitFilt) = -y(hitFilt);
        z(hitFilt) = -z(hitFilt);
        for ee=1:length(stackScanRange)-1
            if ismember(ee,hitFilt)
                color = 'red'; % Between those two time steps the filter wall is penetrated
            else
                color = 'blue'; % Everyting ok here.
            end
            plot3(x(ee:ee+1), y(ee:ee+1), z(ee:ee+1),color)
        end
        % Put a marker (point) at the last vector position to indicate the
        % particle direction.
        plot3(x(length(stackScanRange)),...
            y(length(stackScanRange)),...
            z(length(stackScanRange)),[color '.'],'MarkerSize',5)
    else % everything ok
        j = j+1;
        color = 'blue';
        plot3(x,y,z,color)
        % Put a marker (point) at the last vector position to indicate the
        % particle direction.
        plot3(x(length(stackScanRange)),...
            y(length(stackScanRange)),...
            z(length(stackScanRange)),[color '.'],'MarkerSize',5)
    end
end
grid on
view(-30,45)
