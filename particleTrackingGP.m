%%% changelog.


%% User input
ispc= 0; %0 if linux, 1 if winOS
% linux work:
pathBase = '/home/gp/Pesch/Ph.D./unshared/APS_Tomographie/Exp';
pathRest = '/Thres=0.0005,MedRad=2/';
% mac home:
% stackBasePath = '/Volumes/Untitled/APS_Tomography/Exp041/Thres=0.0005,MedRad=0/bin_part_proj_0121_scan_'
% win OS:
%pathBase = 'f:\aps_2-bm_2016_08\exp';
% pathRest='_d07t05n03_v0.068mlmin-1_c02.5e-4w%_multiplescantest_dimax_2x_100mm_0.4msecexptime_360degpersec_rolling_100umluag_1mmc2mmglass_2.657mrad_usarm1.1502_monoy_-10.0_ahutch\';
% ATTENTION: When changing experimentID also change the beginning of pathRest!
experimentID=41;
stackScanRange = 1:5; % actual number of time steps

% specify the search radius 
r0 = 5; % 15 [px]
dr = 2.5; % 5 [px]
rmax = 20; % 50 [px] equals (rmax-r0)/dr = 7 iterations
cone = 1; % cone search or regular search?
boolPos = 1; %do particles move in positive z direction (1) or negative z direction (0)?

% particle volume, comes from particle size distribution plots and known
% particle diameter (d50 = 18�m = 3.6px)
minParticleSize = 7; %[px^3] equals d = 2.5px
maxParticleSize = 100; %[px^3] equals d = 5.8px

%% Calculation
% iteration through TIME STEPS
for i=stackScanRange
    disp(['Time step is ' num2str(find(i==stackScanRange))]);
    
    % create stack base file path for real scans
    % stackBaseFile = pathCreation(pathBase, pathRest, experimentID, i, ispc);
    
    % create stack base file path for benchmark images
    stackBaseFile = ['benchmarkPics/time' num2str(i) '/im*.tif'];
    stack = readtimeseries(stackBaseFile,'',[],0,0);

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
    
    % first time step t=t_1?
    if i == stackScanRange(1)
        % initialize particle array
        % three-dimensional; dim 1 = [x y z], dim 2 = particle number ii,
        % dim 3 = time step i
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
        %                   at least two particles from ftime step i-1 
        %                   found the same particle in time step i
        %                   o o ---> o  (bad)
        % Particles that do not appear in one of those three arrays had no
        % match in time step i. The radius is increased to deal with those.
        
        found = int16.empty;
        skip = int16.empty;
        multipleFound = int16.empty;
            
        particlesBuffer=zeros(4,size(particles,2));
        
        for ri=r0:dr:rmax
            disp(['Search radius ri =   ' num2str(ri) ]);
            disp(['Considering          ' num2str(length(nonzeros(particles(1,:,i-1)))) ' / ' num2str(size(particles,2)) ' particles' ])
            disp(['Possible partners    ' num2str(length(msr.Id)) ]);
            
            for ii=1:size(particles,2) % iterate through all particles found in t=t_1
                % skip all particles that already have partners (1 or more)
                if ismember(ii,horzcat(found,skip,multipleFound))
                    continue
                end
                
                % skip all particles that didn't have a partner the
                % timestep before
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
                
                [nearbyIds,nearbyCenters] = findNearbyParticles(particles(:,ii,i-1),msrCenters,ri,cone,boolPos);
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
            % clarify which particles from time step i are found.
            for kk=1:max(particlesBuffer(1,:)) % max nearbyIds in particlesBuffer
                % The position of particlesBuffer equals ii (here tmp), the
                % value equals nearbyIds (here kk).
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
            % Sort found t.
%             found = found';
%             found = sortrows(found,1);
%             found = found';
            found = sort(found);
            
            
            % Remove it from measure array by putting zeros in the
            % array at the specific row.
            [~,msrCenters]=treatArray(msrCenters,'remove',particlesBuffer(1,found));
            
            disp(['Found partners in ri ' num2str(size(found,2)) ]);
            disp(['Too many partners    ' num2str(size(skip,2)) ]);
            disp(['Found same particle  ' num2str(size(multipleFound,2)) ]);
            disp(['No match found       ' num2str(length(nonzeros(particles(1,:,i-1)))-size(multipleFound,2)-size(skip,2)-size(found,2)) ]);
            disp('---');
        end
        % Save new particle position.
        particles(:,found,i) = particlesBuffer(2:4,found);
        disp('=====================================================');
    end     
end

%% Display
% strip array of all particles that contain 0 coordinates:
figure;
hold on
j = 1;
for aa=1:size(particles,2)
    % rearrange
    x = reshape(particles(1,aa,:),[length(stackScanRange) 1]);
    y = reshape(particles(2,aa,:),[length(stackScanRange) 1]);
    z = reshape(particles(3,aa,:),[length(stackScanRange) 1]);
    if any(x + y + z == 0)
        continue
    else
        j= j+1;
        plot3(x,y,z)
        % Put a marker (point) at the last vector position to indicate the
        % particle direction.
        plot3(x(length(stackScanRange)), ...
            y(length(stackScanRange)), ...
            z(length(stackScanRange)),['r.'],'MarkerSize',5)
    end
    %if j>=50
    %    break
    %end
end
