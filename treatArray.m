function [msrIDs, msrCenters ] = treatArray(array, str, id)
% This function alters the msr array in the way specified by string 'str'.
% The integer value id is used to replace row no. id in the msr array with zeros.
    if(strcmp(str,'sort'))
        % array = msr
        msrCenters = zeros(4,size(array,1));
        msrCenters(1,:) = array.id;
        msrCenters(2:4,:) = array.center;
        msrCenters = msrCenters';   % Transpose msrCenters to use sortrows().
        msrCenters = sortrows(msrCenters,1);
        % The IDs are stored in this vector in ascending order. They 
        % correspond to msrCenters
        msrIDs = msrCenters(:,1)';
        % Now msrCenters has 3 rows and as many columns as there are IDs.
        msrCenters=msrCenters(:,2:4)';
    elseif(strcmp(str,'remove'))
        % array = msrCenters, id = nearbyIds
        array(:,id)=0;
        msrCenters=array;
        msrIDs=[];  % not used
    end
end

