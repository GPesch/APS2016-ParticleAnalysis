function [ids, centers] = findNearbyParticles(pos, msrCenters, r, cone, boolPos)
    % check wether simple or with cone. Simple just looks for particles in
    % spheres around the original position.
    % cone looks in smaller spheres but also in a cone which is open
    % towards negative z (thus in flow direction!)
    % boolPos decides whether particles flow in positive z direciton or
    % negative z direction

    % Look for other particles within circle (radius r) around particle ii
    % (specified by its postion pos).
    
    
    ids = inCircle(pos,msrCenters,r);
    %ids = [];
    centers=[];
    
    if cone
        idsFromCone = inCone(pos,msrCenters, 3/4*pi, 4*r, boolPos);
        
        % just add the ids that are NOT double.
        % find the index of double ids in Cone array:
        [~, ~, idIdCone] = intersect(ids,idsFromCone);
        % horizontally stack ids from circle and the ids from cone that are
        % not already in the id list:
        idsFromCone(idsFromCone(idIdCone))=[];
        ids = horzcat(ids,idsFromCone);
    end
    
   
    
    if(~isempty(ids))
        centers = msrCenters(:,ids);
    end
    
        
end


function [boolReturnVector] = inCircle(pos, msrCenters, r)
    % returns all INDICES OF THE MEASURE ARRAY (as opposed to real particle
    % ids msr.Id) that are inside the circle

    boolReturnVector = find( (pos(1) - msrCenters(1,:) ).^2 ...
                           + (pos(2) - msrCenters(2,:) ).^2 ...
                           + (pos(3) - msrCenters(3,:) ).^2 ...
                           < r^2 );
end

function [boolReturnVector] = inCone(pos, msrCenters, aperture, h, boolPos)
    % returns all INDICES OF THE MEASURE ARRAY that are inside a cone which
    % is open downwards and whose apex is at the original particle
    % position
    
    
    if boolPos
        zCondition = msrCenters(3,:)-pos(3);
    else
        zCondition = msrCenters(3,:)+pos(3);
    end
        
    boolReturnVector = find( ((pos(1) - msrCenters(1,:)).^2 + ...
                              (pos(2) - msrCenters(2,:)).^2) ...
                              *cos(aperture/2)^2 -...
                              (pos(3) - msrCenters(3,:)).^2 ...
                              *sin(aperture/2)^2 ...
                              < 0 & ...
                              (zCondition) > 0 & ...
                              (zCondition) < h ...
                            );
end