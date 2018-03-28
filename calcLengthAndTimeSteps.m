function [ idsOfInterest, pathLength ] = calcLengthAndTimeSteps( particles, relevTimeSteps, stackScanRange, pathBase, pathRest, experimentID )
    % Check in how many time steps matches for a particle of the first time
    % step were found. Based on that, calculate the path length. Also save
    % particles that are found in e.g. 10 instead of max 20 time steps.
    %
    % Input values: the 'particles' array (dim1 = [x y z] -> negative values 
    % indicate a path that penetrates the solid filter, dim2 = particle
    % ids of first time step, dim3 = time steps); relevant time steps (e.g.
    % you are interested in particles that were matched in 10 time steps out 
    % of max 20); pathBase, pathRest & experimentID to create path of files.
    %
    % Return values: idsOfInterest and pathLength. The actual results are 
    % also written as txt-files in the specified folder.

    % initialize variable
    idsOfInterest = int16.empty;

    for ff = 1:size(particles,2)
        % rearrange
        x = reshape(particles(1,ff,:),[length(stackScanRange) 1]);
        y = reshape(particles(2,ff,:),[length(stackScanRange) 1]);
        z = reshape(particles(3,ff,:),[length(stackScanRange) 1]);

        % check if there is a time step t_i within the relevant time steps
        % (user input) where particles were not found (i.e. x(i) & y(i) & z(i) are all 0).
        % Here [x y z] gives an array with 3 cols and as many rows as relevant
        % time steps. ismember() checks for zeros, the cols are summed up and
        % if there is no sum=3 (again checked by ismember()) it means that all
        % relevant time steps have a position/an entry, i.e. they were matched.
        if ~any(ismember(sum(ismember([x(1:relevTimeSteps) y(1:relevTimeSteps) z(1:relevTimeSteps)],0),2),3))
            idsOfInterest = horzcat(idsOfInterest,ff); % add interesting ID
        end
    end

    % Create array with interesting particles (particles that were found in at
    % least x (-> user input) time steps).
    particlesOfInterest = particles(:,idsOfInterest,:);
    %save
    dlmwrite([pathBase num2str(experimentID,'%03.0f') pathRest 'particlesTracked_in' num2str(relevTimeSteps) 'timeStepsMin_particleTrackingGP_0_0_1_FiltEval0_2_1.txt'],particlesOfInterest,'delimiter','\t','precision',8);

    %calculate path length
    pathLength = zeros(size(idsOfInterest,2),3); % 3 cols
    pathLength(:,1) = idsOfInterest'; % 1st col is filled with idsOfInterest

    for gg = 1:size(idsOfInterest,2)
        % rearrange
        % abs() is necessary in order to deal with possible negative
        % positions (position is negative if path penetrates the filter -> 
        % see main script). abs() makes negative entries positive.
        x = abs(reshape(particles(1,idsOfInterest(gg),:),[length(stackScanRange) 1]));
        y = abs(reshape(particles(2,idsOfInterest(gg),:),[length(stackScanRange) 1]));
        z = abs(reshape(particles(3,idsOfInterest(gg),:),[length(stackScanRange) 1]));

        % 2nd col is filled with actual time steps where particles were matched.
        if isempty(min(find((x + y + z)==0,1))) % no zeros...particle in all time steps found
            pathLength(gg,2) = 20;
        else
            pathLength(gg,2) = min(find((x + y + z)==0,1)) - 1;
        end

        % 3rd col is filled with path length.
        pathLeng = 0;
        for ggg = 2:pathLength(gg,2)
            pathLeng = pathLeng + norm([(x(ggg) - x(ggg - 1)) (y(ggg) - y(ggg - 1)) (z(ggg) - z(ggg - 1))]);
        end
        pathLength(gg,3) = pathLeng;
    end
    %save
    dlmwrite([pathBase num2str(experimentID,'%03.0f') pathRest 'particlesTracked_pathLength_particleTrackingGP_0_0_1_FiltEval0_2_1.txt'],pathLength,'delimiter','\t','precision',8);

end

