[rawim, XYZ2Cam, wbcoeffs] = readdng('RawImage.DNG');

[Csrgb, Clinear, Cxyz, Ccam] = dng2rgb(rawim, XYZ2Cam, wbcoeffs, 'rggb', 'nearest', 4000, 6000);

%figure 1 has all 4 images
figure(1)
subplot(2,2,1);
imshow(Ccam);
title('Ccam');

subplot(2,2,2);
imshow(Cxyz);
title('Cxyz');

subplot(2,2,3);
imshow(Clinear);
title('Clinear');

subplot(2,2,4);
imshow(Csrgb);
title('Csrgb');

%figure 2 is Csrgb
figure(2)
imshow(Csrgb);

%figure 3 has Ccam's histograms
figure(3);
subplot(1,3,1);
imhist(Ccam(:,:,1));
title('Red channel Histogram');
subplot(1,3,2);
imhist(Ccam(:,:,2));
title('Green channel Histogram');
subplot(1,3,3);
imhist(Ccam(:,:,3));
title('Blue channel Histogram');

%figure 4 has Cxyz's histograms
figure(4);
subplot(1,3,1);
imhist(Cxyz(:,:,1));
title('Red channel Histogram');
subplot(1,3,2);
imhist(Cxyz(:,:,2));
title('Green channel Histogram');
subplot(1,3,3);
imhist(Cxyz(:,:,3));
title('Blue channel Histogram');

%figure 5 has Clinear's histograms
figure(5);
subplot(1,3,1);
imhist(Clinear(:,:,1));
title('Red channel Histogram');
subplot(1,3,2);
imhist(Clinear(:,:,2));
title('Green channel Histogram');
subplot(1,3,3);
imhist(Clinear(:,:,3));
title('Blue channel Histogram');

%figure 6 has Csrgb's histograms
figure(6);
subplot(1,3,1);
imhist(Csrgb(:,:,1));
title('Red channel Histogram');
subplot(1,3,2);
imhist(Csrgb(:,:,2));
title('Green channel Histogram');
subplot(1,3,3);
imhist(Csrgb(:,:,3));
title('Blue channel Histogram');
