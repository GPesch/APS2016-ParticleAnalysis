function [ids, centers] = findNearbyParticles(pos, msrCenters, r)
    % check wether simple or with cone. Simple just looks for particles in
    % spheres around the original position.
    % cone looks in smaller spheres but also in a cone which is open
    % towards negative z (thus in flow direction!)

    % Look for other particles within circle (radius r) around particle ii
    % (specified by its postion pos). Returns id(s).
    ids = inCircle(pos,msrCenters,r);
    centers=[];
   
    
    if(~isempty(ids))
        centers = msrCenters(:,ids);
    end
    
        
end


function [boolReturnVector] = inCircle(pos, msrCenters, r)
    % returns all INDICES OF THE MEASURE ARRAY (as opposed to real particle
    % ids msr.Id) that are inside the circle. The measure array olds the
    % positions of the CURRENT time step. The found id will be stored in
    % particlesBuffer in order to check if multiple matches were found. The
    % ids have nothing to do with the loop variable ii (particle id of the
    % FIRST time step) which is the input value for the main function findNearbyParticles().

    % find() returns the position of the match, i.e. the id (because
    % msrCenters is sorted by ids).
    boolReturnVector = find( (pos(1) - msrCenters(1,:) ).^2 ...
                           + (pos(2) - msrCenters(2,:) ).^2 ...
                           + (pos(3) - msrCenters(3,:) ).^2 ...
                           < r^2 );
end