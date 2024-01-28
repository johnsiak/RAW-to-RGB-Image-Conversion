function [Csrgb, Clinear, Cxyz, Ccam] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, ...
    bayertype, method, M, N)

[Mo, No] = size(rawim);

pattern = zeros(2,2); 

mask = wbcoeffs(2)*ones(Mo,No); %Initialize to all green values

switch bayertype
    case 'rggb'
        mask(1:2:end, 1:2:end) = wbcoeffs(1); %r
        mask(2:2:end, 2:2:end) = wbcoeffs(3); %b
        pattern = [1 2; 2 3]; %1 for red, 2 for green, 3 for blue
    case 'bggr'
        mask(2:2:end, 2:2:end) = wbcoeffs(1); %r
        mask(1:2:end, 1:2:end) = wbcoeffs(3); %b
        pattern = [3 2; 2 1]; %1 for red, 2 for green, 3 for blue
    case 'grbg'
        mask(1:2:end, 2:2:end) = wbcoeffs(1); %r
        mask(2:2:end, 1:2:end) = wbcoeffs(3); %b
        pattern = [2 1; 3 2]; %1 for red, 2 for green, 3 for blue
    case 'gbrg'
        mask(2:2:end, 1:2:end) = wbcoeffs(1); %r
        mask(1:2:end, 2:2:end) = wbcoeffs(3); %b
        pattern = [2 3; 1 2]; %1 for red, 2 for green, 3 for blue
end

%white balance
balanced_rawim = rawim .* mask;

rgbim = zeros(Mo,No,3);

if(method == "nearest")
    %near neighbor method
    switch bayertype
        case 'rggb'
            %we calculate the RGB values of each pixel in a 2x2 block
            for j = 1:2:Mo-1
                for i = 1:2:No-1
                    %red
                    rgbim(j,i,1) = balanced_rawim(j,i);
                    rgbim(j,i+1,1) = balanced_rawim(j,i);
                    rgbim(j+1,i,1) = balanced_rawim(j,i);
                    rgbim(j+1,i+1,1) = balanced_rawim(j,i);
                    %green 1
                    rgbim(j,i,2) = balanced_rawim(j,i+1); 
                    rgbim(j,i+1,2) = balanced_rawim(j,i+1);
                    %green 2
                    rgbim(j+1,i,2) = balanced_rawim(j+1,i);
                    rgbim(j+1,i+1,2) = balanced_rawim(j+1,i);
                    %blue
                    rgbim(j,i,3) = balanced_rawim(j+1,i+1);
                    rgbim(j,i+1,3) = balanced_rawim(j+1,i+1);
                    rgbim(j+1,i,3) = balanced_rawim(j+1,i+1);
                    rgbim(j+1,i+1,3) = balanced_rawim(j+1,i+1);
                end
            end
        case 'bggr'
            for j = 1:2:Mo-1
                for i = 1:2:No-1
                    %red
                    rgbim(j,i,1) = balanced_rawim(j+1,i+1);
                    rgbim(j,i+1,1) = balanced_rawim(j+1,i+1);
                    rgbim(j+1,i,1) = balanced_rawim(j+1,i+1);
                    rgbim(j+1,i+1,1) = balanced_rawim(j+1,i+1);
                    %green 1
                    rgbim(j,i,2) = balanced_rawim(j,i+1); 
                    rgbim(j,i+1,2) = balanced_rawim(j,i+1);
                    %green 2
                    rgbim(j+1,i,2) = balanced_rawim(j+1,i);
                    rgbim(j+1,i+1,2) = balanced_rawim(j+1,i);
                    %blue
                    rgbim(j,i,3) = balanced_rawim(j,i);
                    rgbim(j,i+1,3) = balanced_rawim(j,i);
                    rgbim(j+1,i,3) = balanced_rawim(j,i);
                    rgbim(j+1,i+1,3) = balanced_rawim(j,i);
                end
            end
        case 'grbg'
            for j = 1:2:Mo-1
                for i = 1:2:No-1
                    %red
                    rgbim(j,i,1) = balanced_rawim(j,i+1);
                    rgbim(j,i+1,1) = balanced_rawim(j,i+1);
                    rgbim(j+1,i,1) = balanced_rawim(j,i+1);
                    rgbim(j+1,i+1,1) = balanced_rawim(j,i+1);
                    %green 1
                    rgbim(j,i,2) = balanced_rawim(j,i); 
                    rgbim(j,i+1,2) = balanced_rawim(j,i);
                    %green 2
                    rgbim(j+1,i,2) = balanced_rawim(j+1,i+1);
                    rgbim(j+1,i+1,2) = balanced_rawim(j+1,i+1);
                    %blue
                    rgbim(j,i,3) = balanced_rawim(j+1,i);
                    rgbim(j,i+1,3) = balanced_rawim(j+1,i);
                    rgbim(j+1,i,3) = balanced_rawim(j+1,i);
                    rgbim(j+1,i+1,3) = balanced_rawim(j+1,i);
                end
            end
        case 'gbrg'
            for j = 1:2:Mo-1
                for i = 1:2:No-1
                    %red
                    rgbim(j,i,1) = balanced_rawim(j+1,i);
                    rgbim(j,i+1,1) = balanced_rawim(j+1,i);
                    rgbim(j+1,i,1) = balanced_rawim(j+1,i);
                    rgbim(j+1,i+1,1) = balanced_rawim(j+1,i);
                    %green 1
                    rgbim(j,i,2) = balanced_rawim(j,i); 
                    rgbim(j,i+1,2) = balanced_rawim(j,i);
                    %green 2
                    rgbim(j+1,i,2) = balanced_rawim(j+1,i+1);
                    rgbim(j+1,i+1,2) = balanced_rawim(j+1,i+1);
                    %blue
                    rgbim(j,i,3) = balanced_rawim(j,i+1);
                    rgbim(j,i+1,3) = balanced_rawim(j,i+1);
                    rgbim(j+1,i,3) = balanced_rawim(j,i+1);
                    rgbim(j+1,i+1,3) = balanced_rawim(j,i+1);
                end
            end  
    end 

elseif(method == "bilinear")
    %bilinear method

    %each pixel corresponds to a colour 
    %1 for red, 2 for green, 3 for blue
    fullpattern = zeros(Mo,No);
    
    %colour assignment to every pixel
    for j = 1:2:Mo-1
        for i = 1:2:No-1
            fullpattern(j:j+1,i:i+1) = pattern;
        end
    end
   
    %we extend our matrices because we want to apply the bilinear method
    %to all pixels including the border.
    %by using this approach we ensure that we do not exceed
    %the boundaries of the matrix
    fullpattern_extended = zeros(Mo+2,No+2);
    fullpattern_extended(2:Mo+1,2:No+1) = fullpattern;

    balanced_rawim_extended = zeros(Mo+2,No+2);
    balanced_rawim_extended(2:Mo+1,2:No+1) = balanced_rawim;

    for j = 2:Mo
        for i = 2:No
            value1 = 0;
            count1 = 0;
            value2 = 0;
            count2 = 0;

            %neighbors is a 2x8 array with information about the 8 neighbors
            %the first row is the colour of the pixel
            %the second row is the rawim value of the pixel
            neighbors = [fullpattern_extended(j-1,i-1), fullpattern_extended(j-1,i), fullpattern_extended(j-1,i+1), ...
                fullpattern_extended(j,i-1), fullpattern_extended(j,i+1), ...
                fullpattern_extended(j+1,i-1), fullpattern_extended(j+1,i), fullpattern_extended(j+1,i+1);

                balanced_rawim_extended(j-1,i-1), balanced_rawim_extended(j-1,i), balanced_rawim_extended(j-1,i+1), ...
                balanced_rawim_extended(j,i-1), balanced_rawim_extended(j,i+1), ...
                balanced_rawim_extended(j+1,i-1), balanced_rawim_extended(j+1,i), balanced_rawim_extended(j+1,i+1)];

            if fullpattern_extended(j,i) == 1 %if pixel is red
                for k = 1:8 % 8 neighbors
                    if(neighbors(1,k) == 2) %if neighbor is green
                        value1 = value1 + neighbors(2,k);
                        count1 = count1 + 1;
                    elseif(neighbors(1,k) == 3) %if neighbor is blue
                        value2 = value2 + neighbors(2,k);
                        count2 = count2 + 1;
                    end
                end

                rgbim(j-1,i-1,1) = balanced_rawim_extended(j,i); %red
                rgbim(j-1,i-1,2) = value1/count1; %green
                rgbim(j-1,i-1,3) = value2/count2; %blue

            elseif fullpattern_extended(j,i) == 2 %if pixel is green
                for k = 1:8 % 8 neighbors
                    if(neighbors(1,k) == 1) %if neighbor is red
                        value1 = value1 + neighbors(2,k);
                        count1 = count1 + 1;
                    elseif(neighbors(1,k) == 3) %if neighbor is blue
                        value2 = value2 + neighbors(2,k);
                        count2 = count2 + 1;
                    end
                end

                rgbim(j-1,i-1,1) = value1/count1; %red
                rgbim(j-1,i-1,2) = balanced_rawim_extended(j,i); %green
                rgbim(j-1,i-1,3) = value2/count2; %blue

            elseif fullpattern_extended(j,i) == 3 %if pixel is blue
                for k = 1:8 % 8 neighbors
                    if(neighbors(1,k) == 1) %if neighbor is red
                        value1 = value1 + neighbors(2,k);
                        count1 = count1 + 1;
                    elseif(neighbors(1,k) == 2) %if neighbor is green
                        value2 = value2 + neighbors(2,k);
                        count2 = count2 + 1;
                    end
                end

                rgbim(j-1,i-1,1) = value1/count1; %red
                rgbim(j-1,i-1,2) = value2/count2; %green
                rgbim(j-1,i-1,3) = balanced_rawim_extended(j,i); %blue

            end            
        end
    end
end

%sampling factors
row_ratio = (Mo-1)/(M-1);
col_ratio = (No-1)/(N-1);

rgbim_resized = zeros(M, N, 3);

%we assign the 4 corner-pixels of the rgbim to rgbim_resized
rgbim_resized(1, 1, :) = rgbim(1, 1, :);
rgbim_resized(1, end, :) = rgbim(1, No, :);
rgbim_resized(end, 1, :) = rgbim(Mo, 1, :);
rgbim_resized(end, end, :) = rgbim(Mo, No, :);

%sampling
for j = 2:M-1
    for i = 2:N-1
        x = round((j-1)*row_ratio) + 1;
        y = round((i-1)*col_ratio) + 1;
        rgbim_resized(j,i,:) = rgbim(x,y,:);
    end
end

Ccam = rgbim_resized;

Cxyz = zeros(M,N,3);

for j = 1:M
    for i = 1:N
        tmp = reshape(Ccam(j,i,:), [3 1]); %reshape to 2D array
        Cxyz(j,i,:) = XYZ2Cam\tmp; %XYZ2Cam^(-1) * tmp
    end
end

XYZ2RGB = [3.2406 -1.5372 -0.4986;
           -0.9689 1.8758 0.0415; 
            0.0557 -0.2040 1.0570];

%normalizing every row so that each row's sum is 1
XYZ2RGB = XYZ2RGB./repmat(sum(XYZ2RGB,2),1,3); 

Clinear = zeros(M,N,3);

for j = 1:M
    for i = 1:N
        tmp = reshape(Cxyz(j,i,:), [3 1]); %reshape to 2D array
        Clinear(j,i,:) = XYZ2RGB*tmp;
    end
end

%brightness correction
grayim = rgb2gray(Clinear);
grayscale = 0.25/mean(grayim(:));
bright_srgb = min(1,Clinear*grayscale);

%gamma correction
Csrgb = bright_srgb.^(1/2.2);

end