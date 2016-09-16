function [ids, centers] = findNearbyParticles(pos, msrCenters, r, cone)
    % check wether simple or with cone. Simple just looks for particles in
    % spheres around the original position.
    % cone looks in smaller spheres but also in a cone which is open
    % towards negative z (thus in flow direction!)

    % Look for other particles within circle (radius r) around particle ii
    % (specified by its postion pos).
    ids = inCircle(pos,msrCenters,r)
    centers=[];
   
    
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