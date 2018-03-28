function [ ] = drawPathInImage( particles, stackScanRange, idsOfInterest, imageZ, pathBase, pathRest, experimentID )
    % Draws the paths of found particles in a 3D image and adds spheres to
    % the last coordinates in order to indicate the direction.

    particlesOfInterest = particles(:,idsOfInterest,:);
    pathImage = newim([2016,2016,imageZ], 'uint8');
    
    startPoint = uint16.empty;
    endPoint = uint16.empty;
    sphereCoordinates = uint16.empty;
    sphereRadius = 3; % [px]
        
    % loop through particlesOfInterest
    for hh = 1:size(idsOfInterest,2)
        % rearrange
        x = reshape(particlesOfInterest(1,hh,:),[length(stackScanRange) 1]);
        y = reshape(particlesOfInterest(2,hh,:),[length(stackScanRange) 1]);
        z = reshape(particlesOfInterest(3,hh,:),[length(stackScanRange) 1]);
        
        % loop through time steps
        for hhh = 1:(length(stackScanRange) - 1)
            % [x y z] might be filled with '0's at some point (if no
            % particles have been matched after that time). In order not
            % to add the coordinate [0 0 0] to startPoint and endPoint (a
            % line would be drawn to that point) break here.
            if (round(x(hhh + 1) + y(hhh + 1) + z(hhh + 1)) == 0)
                % The startpoint (hhh) might not be zero but if the 
                % endpoint (hhh + 1) is zero save the last position of the 
                % path, then break and handle the next particle. 
                sphereCoordinates = vertcat(sphereCoordinates, round( [x(hhh) y(hhh) z(hhh)] ));
                break
            end
            
            % If the last time step was reached -> save position for sphere
            if hhh == (length(stackScanRange) - 1)
                sphereCoordinates = vertcat(sphereCoordinates, round( [x(hhh) y(hhh) z(hhh)] ));
            end
            
            startPoint = vertcat(startPoint, round( [x(hhh) y(hhh) z(hhh)] ));
            endPoint = vertcat(endPoint, round( [x(hhh + 1) y(hhh + 1) z(hhh + 1)] ));
            % They have the form of (?? rows, 3 cols = [x1 y1 z1; x2 y2 z2;...] ):
            %
            % x1 y1 z1
            % x2 y2 z2
            % .  .  .
            % .  .  .
            % .  .  .
        end
    end
    
    intensity = 255;
    % First make spheres, then draw lines. Order is irrelevant though.
    pathImage = addSpheres(pathImage, sphereCoordinates, sphereRadius, intensity);
    % last parameter = intensity -> uint8 -> max intensity = 255.
    % startPoint and endPoint must be integers (hence round()).
    % the start point of line 2 can be different from the end point of line 1 
    % (i.e. several not connected lines can be drawn).
    pathImage = drawline(pathImage, startPoint, endPoint, intensity);
    % Finally, dilate objects (paths are only one px thick).
    iterations = 3;
    pathImage = dip_image(uint8((bdilation(pathImage>0, iterations, -3, 0) * 255)));
    
    % Save image as ics-file (can be converted to tiff-stack using ImageJ).
    writeim(pathImage,[pathBase num2str(experimentID,'%03.0f') pathRest 'particlePath.ics'],'ICSv2',0,[]);
end

function [ image ] = addSpheres(image, coordinates, radius, intensity)
    % Add spheres to end points of particle paths to indicate direction of 
    % particle movement.
    
    % Create array (matrix) of image -> saves A LOT of computational time (-> for-loops are critical).
    % ATTENTION: imageAr(x+1,y+1,z+1) == image(y,x,z)
    %   -> x and y are switched
    %   -> image starts with index 0, array with index 1
    imageAr = dip_array(image);
    % loop through particles
    for pp = 1:size(coordinates,1)
        % Avoid looping through whole image/array (time intensive). Hence loop
        % through subimage/subarray of the size of the sphere.
        x1 = coordinates(pp,1) - radius;
        x2 = coordinates(pp,1) + radius;
        y1 = coordinates(pp,2) - radius;
        y2 = coordinates(pp,2) + radius;
        z1 = coordinates(pp,3) - radius;
        z2 = coordinates(pp,3) + radius;
        
        % Handle problematic coordinates (e.g. z1 < 0)
        if x1 < 0 
            x1 = 0;
        end
        if y1 < 0
            y1 = 0;
        end
        if z1 < 0
            z1 = 0;
        end
        if x2 > (size(image,1) - 1)
            x2 = size(image,1) - 1;
        end
        if y2 > (size(image,2) - 1)
            y2 = size(image,2) - 1;
        end
        if z2 > (size(image,3) - 1)
            z2 = size(image,3) - 1;
        end
        
        % loop through subimage/subarray in order to draw sphere.
        for x = x1:x2
            for y = y1:y2
                for z = z1:z2
                    % draw sphere
                    if ((x - coordinates(pp,1))^2 + (y - coordinates(pp,2))^2 + (z - coordinates(pp,3))^2 <= radius^2)
                        imageAr(y + 1, x + 1, z + 1) = intensity;
                        % if image:
                        % image(x,y,z) = intensity;
                    end
                end
            end
        end
        % Done with particle pp.
    end
    % Done with all particles.
    
    % Change data type back to dip_image.
    image = dip_image(imageAr);
end

