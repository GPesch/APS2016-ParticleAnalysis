function [imageStack] = placeParticle(imageStack,pos,dPart)
    % considering the following range of pixels in image stack
    minMax = horzcat(floor(pos-dPart),ceil(pos+dPart));
    for ix = minMax(1,1):minMax(1,2)
        for iy = minMax(2,1):minMax(2,2)
            for iz = minMax(3,1):minMax(3,2)
                if pixelDistSqr([ix-0.5; iy-0.5; iz-0.5],pos) < (dPart/2)^2
                    imageStack(ix,iy,iz)=1;
                end
            end
        end
    end
end

function [distance] = pixelDistSqr(pos1,pos2)
    distance = (pos1(1)-pos2(1))^2+(pos1(2)-pos2(2))^2+(pos1(3)-pos2(3))^2;
end