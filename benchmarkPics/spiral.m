% user input

pictSize = [500 500 200] %x y z
noParticles = 30
dPart = 3 % px
spiralRadMax = 25 % px % maximum value
t = 0:1:10
spiralPeriod = 10
velocity = 6


% initialize start point
possibleStart(1,:) = [0+spiralRadMax+dPart/2 pictSize(1)-spiralRadMax-dPart/2]
possibleStart(2,:) = [0+spiralRadMax+dPart/2 pictSize(2)-spiralRadMax-dPart/2]
possibleStart(3,:) = [dPart/2 pictSize(3)/2]

% initialize particle vector

particles = zeros(3,noParticles,length(t))

for i=1:noParticles
    spiralRad = spiralRadMax*rand();
    particles(:,i,1) = possibleStart(:,1) + (possibleStart(:,2) - possibleStart(:,1)).*rand(3,1); %x
    for ii=2:length(t)
        particles(:,i,ii) = [spiralRad*cos(2*pi*t(ii)/spiralPeriod*velocity); spiralRad*sin(2*pi*t(ii)/spiralPeriod*velocity); t(ii)*velocity] + particles(:,i,1) + [spiralRad;0;0];
    end
end

for i=1:size(particles,3)
    timeFolder = ['time' num2str(i)]
    mkdir(timeFolder)
    
    % create the empty iamge matrix
    imageStack = zeros(pictSize);
    ti = t(i);
    for ii=1:size(particles,2)
        imageStack = placeParticle(imageStack,particles(:,ii,i),dPart);
    end
    
    % iterate stack and write image
    for ii = 1:size(imageStack,3)
        fileName = [timeFolder '/im' num2str(ii) '.tif']
        imwrite(imageStack(:,:,ii), fileName, 'tiff')
    end
end