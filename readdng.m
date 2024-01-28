function [rawim, XYZ2Cam, wbcoeffs] = readdng(filename)
    
%open .dng
obj = Tiff(filename, 'r');
offsets = getTag(obj, 'SubIFD');
setSubDirectory(obj, offsets(1));
rawim = read(obj);
close(obj);
    
%meta-data
meta_info = imfinfo(filename);

%width and height of the image (the useful part of array rawim)
width = meta_info.SubIFDs{1}.DefaultCropSize(1);
height = meta_info.SubIFDs{1}.DefaultCropSize(2);

%(x_origin ,y_origin) is the uper left corner of the useful part of the
%sensor and consequently of the array rawim
y_origin = meta_info.SubIFDs{1}.ActiveArea(1) + 1;
x_origin = meta_info.SubIFDs{1}.ActiveArea(2) + 1;

rawim = double(rawim(y_origin:y_origin+height-1, x_origin:x_origin+width-1));

blacklevel = meta_info.SubIFDs{1}.BlackLevel(1);
whitelevel = meta_info.SubIFDs{1}.WhiteLevel;

wbcoeffs = (meta_info.AsShotNeutral).^-1;
wbcoeffs = wbcoeffs/wbcoeffs(2);

XYZ2Cam = meta_info.ColorMatrix2;
XYZ2Cam = reshape(XYZ2Cam, 3, 3)';

%normalizing every row so that each row's sum is 1
%sum(XYZ2Cam, 2) returns a column vector containing the sum of each row
XYZ2Cam = XYZ2Cam./repmat(sum(XYZ2Cam,2),1,3);

%linearizing so that blacklevel is 0 and whitelevel is 1
rawim = (rawim - blacklevel)/(whitelevel - blacklevel);
rawim = max(0,min(rawim,1));

end