function [pathReturn] = pathCreation(pathBase, pathRest, experimentID, scanID, ispc, binWhat)
    % Creates string to load image stack.
    experiment = num2str(experimentID,'%03.0f');
    % Project ID can be calculated from experiment ID. However, after 
    % experiment 73 project ID counts from 0 again.
    if(experimentID<74)
        projectID = experimentID*3-2;
    else
        projectID = (experimentID-74)*3+2;
    end
    project = num2str(projectID,'%04.0f');
    scan = num2str(scanID,'%02.0f');
    
    % care for OS -> '/' or '\' for subfolders
    if ispc
        pathReturn = [pathBase experiment pathRest 'bin_' binWhat '_proj_' project '_scan_' scan '\bin_' binWhat '_proj_' project '_scan_' scan '_00000.tif'];
    else
        pathReturn = [pathBase experiment pathRest 'bin_' binWhat '_proj_' project '_scan_' scan '/bin_' binWhat '_proj_' project '_scan_' scan '_00000.tif'];
    end
end

