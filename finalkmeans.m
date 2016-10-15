%function [inpput_image,nroi_image,roi_image,fuse_image,dc_roi,dc_nroi,D]=seg_kmeans(filename,p)
%%
% close all figure windows, plots
close all
% clear all variables in the current workspace
 clear all;
% reads image data from the compliant DICOM(4,6,9) 
p = 1
filename='i4.dcm';
if(p==1)
    im = dicomread(filename);
    
    % normalize the dicom image to 0-255 range for fast processing
    max_val = max(max(im));
    z = max_val/255;
    im = im./z;
else
%     read jpeg image
    im=imread(filename);
    I=im;
    im=rgb2gray(im);
    im=imresize(im,[512 512]);
end

inpput_image=im;
org_img=im;

I1=im;
im=im2double(im);
 subplot(2,2,1);
% % % %the second parameter controls the display range of a grayscale image
imshow(im,[]), title('input image');
j=im;
%%
%mask size 5 x 5
N = 5;
%to retain original size of image even after filtering
im_pad = padarray(j,[floor(N/2) floor(N/2)]);
%sort elements in ascending or descending order to get median
%convert each sliding NxN block of A into a column of B, with no zero padding
%NOTE im_col is not a column matrix it is a rearranged matrix with 2D as in im_pad
im_col = im2col(im_pad,[N N],'sliding');
%sorts each column of im_pad in ascending order still retains 2D
sorted_cols = sort(im_col,1);
%finds median and returns a row martix containing medians for all columns of the input
%B(x,y) returns the element in B at location (x,y); if y is not specified
%it returns the element at position x in column 1; : implies all columns
med_vector = sorted_cols(floor(N*N/2)+1,:);
%convert result to matrix form to get the MEDIAN FILTERED image
out1 = col2im(med_vector,[N N],size(im_pad),'sliding');
 subplot(2,2,2);
 imshow(out1,[]), title('median filtered image');
%%
% Generate Gaussian mask/kernel
sigma = 2; % Define sigma here; small detects fine features
% one method to create a matrix; ind is a row matrix [-2 -1 0 1 2]
ind = -floor(N/2) : floor(N/2);
% ind=exp(ind);
% create X by repeating ind rowwise and Y by repeating ind columnwise
[X Y] = meshgrid(ind, ind);
% implement Gaussian HP filter formula; h is 5x5
h = exp(-(X.^2 + Y.^2) / (2*sigma*sigma));
% GHPF=1-GLPF
h=1-h;

% Convert filter into a column vector
h = h(:);

% Filter our image
I=out1;
% this conversion is needed because inbuilt command requires compatibility with its arguement; 
% bsxfun requires both C and h to be of the same datatype double
% (convert to higher datatype to prevent data loss) 
I = im2double(I);
I_pad = padarray(I, [floor(N/2) floor(N/2)]);
C = im2col(I_pad, [N N], 'sliding');
% Binary Singleton Expansion Function C = bsxfun(FUNC,A,B) applies the 
% element-by-element binary operation specified by the function to A and B
C_filter = sum(bsxfun(@times, C, h), 1);
% rearranges the row vector B into a matrix of size (MM-M+1)-by-(NN-N+1).
% B must be a vector of size 1-by-(MM-M+1)*(NN-N+1).
out2 = col2im(C_filter, [N N], size(I_pad), 'sliding');
 subplot(2,2,3);
 imshow(out2,[]), title('Gaussian filtered image');
 %%
 workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 24;

%im = dicomread('i9.dcm');
 % normalize the dicom image to 0-255 range for fast processing
    max_val = max(max(im));
    z = max_val/255;
    im = im./z;
% Get the dimensions of the image.
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(im);

% Display the original gray scale image.
figure(5)
subplot(2, 3, 1);
imshow(im,[]);

title('Original Grayscale Image');
% Crop image to get rid of light box surrounding the image
im = im(3:end-3, 4:end-4);
% Threshold to create a binary image
binaryImage = im > 20;
% Get rid of small specks of noise
binaryImage = bwareaopen(binaryImage, 10);
% Display the original gray scale image.
subplot(2, 3, 3);
imshow(binaryImage,[]);

title('Binary Image');
% Seal off the bottom of the head - make the last row white.
binaryImage(end,:) = true;
% Fill the image
binaryImage = imfill(binaryImage, 'holes');
subplot(2, 3, 4);
imshow(binaryImage,[]);

title('Cleaned Binary Image');
% Erode away 15 layers of pixels.
se = strel('disk', 15, 0);
binaryImage = imerode(binaryImage, se);
subplot(2, 3, 5);
imshow(binaryImage, []);

title('Eroded Binary Image');
% Mask the gray image
finalImage = im; % Initialize.
finalImage(~binaryImage) = 0;
figure(1)
 subplot(2,2,3);
 
imshow(finalImage, []);

title('Skull stripped Image');
%%
%//K MEANS ALGORITHM
%out2 is double with max value 0.03... out3 is double with max value 255
%so to calculate histogram more accurately in a time efficient manner we
%convert to uint8 with gray level range 0-255
  out2=finalImage;

max_val=max(max(out2));
out3=out2.*(255/max_val);
im1=uint8(out3);

%choose k according to the number of clusters desired in the output image
k=4;

%histogram calculation
%y axis giving number of pixels of a given gray level(x axis value)
img_hist = zeros(256,1);
%gray level value along the x axis
hist_value = zeros(256,1);

%initialize histogram y axis: total no. of pixels for each gray level
for i=1:256
    img_hist(i)=sum(sum(im1==(i-1)));
    %== performs element by element comparison for equality and returns a
    %vector with same dimensions as the matrices being compared with 1 in
    %place of elements for which the inequality is true and 0 for false
end;
%initialize histogram x axis gray levels 0....255
for i=1:256
    hist_value(i)=i-1;
end;

%cluster initialization
%column matrix to hold the centroids - averege intensity value
cluster = zeros(k,1);
%column matrix to store number of pixels in each cluster 
cluster_count = zeros(k,1);
%select random initial seeds or centroids for the k clusters
for i=1:k
    cluster(i)=uint8(rand*205);
end;
%cluster(i)=randi(255);

%to store the immediately previous centroids so as to check
old = zeros(k,1);

%termination condition
%terminate when the error percentage is within acceptable limits
%i.e. when the sum of differences of two successive iterations is <=k
%k=number of clusters/number of intensity levels in final image
%keep iterating till value obtained from 2 successive iterations is almost
%the same
while (sum(sum(abs(old-cluster))) >k)
    old = cluster;
    %matrix to hold the k-value of the cluster which is closest in
    %intensity to the gray value indicated by (index -1)
    closest_cluster = zeros(256,1);
    %subtracts the constant value cluster(1) from all elements of
    %hist_value which represents (x axis)gray levels from 0 to 255
    %min_distance is a vector that store the intensity difference 
    %initialize 
    min_distance = abs(hist_value-cluster(1));
    
    %calculate the minimum distance to a cluster
    %overall (from any cluster) minimum intensity difference matrix  
    for i=2:k
        min_distance =min(min_distance,  abs(hist_value-cluster(i)));
    end;
    
    %calculate the closest cluster by comparing to check if that is the
	%minimum among all clusters
    %closest_cluster is a column matrix whose entries are the cluster 
    %numbers (k) and indices point to the intensity values that are to
    %belong to that cluster
    %find array elements that meet a certain condition and replace with
    %value specified on RHS
    
    %fill closest_cluster matrix with the k value of the cluster
    %closest in intensity to a given gray level
    for i=1:k
        closest_cluster(min_distance==(abs(hist_value-cluster(i)))) = i;
    end;
    
    %calculate the cluster count
    %calculate the number of pixels belonging to each cluster
    
    %closest_cluster contains 1...20 as array elements indicating the
    %clusters to which the gray levels represented by their indices belong
    %img_hist contain the number of pixels of a particular gray level
    %find an element in the matrix equal to the constant i
    for i=1:k
        cluster_count(i) = sum(img_hist .*(closest_cluster==i));
    end;
    
%calculate the new pixel intensity centroids for all clusters
    for i=1:k
        if (cluster_count(i) == 0)
            %cluster has no pixels at all
            cluster(i) = uint8(rand*255);
        else
            %closest_cluster gives the gray level values belonging to a
            %given cluster i
            cluster(i) = uint8(sum(img_hist(closest_cluster==i).*hist_value(closest_cluster==i))/cluster_count(i));
        %mean intensity or gray level = sum(xn1*yn2)/(n1+n2)
        end;
    end;
    
end;
%output matrix of the same size as input image initialized to zeroes
imresult=uint8(zeros(size(im1)));
%cluster matrix has pixel intensity centroids
%closest_cluster contains the k-value of the cluster 1,2...20
%after this operation imresult holds intensity centroid values
%example: all 0s in the image will be replaced by the intensity centroid
%value of the cluster to which it(gray level 0) belongs to
%each gray level value of the input image is replaced by the intensity
%centroid of the cluster to which it belongs, in the output image
for i=1:256
    imresult(im1==(i-1))=cluster(closest_cluster(i));
end;
figure(1)
subplot(2,2,4);
 imshow(imresult,[]);
 title('k means output');
 imtool(imresult,[]);
%%
% obtaining ROI and NROI
im_bin=im2bw(imresult);
[r,c]=find(im_bin==1);
rc=[r c];
simg=im;

im=im2uint8(im);
max_val = max(max(simg));
z = max_val/255;
simg = simg./z;
simg=uint8(simg);
orgimg=simg;

for j=1:(numel(rc)/2)
    simg(r(j),c(j))=0;
end
 figure, subplot(2,1,1); 
 imshow(simg,[]), title('NROI-Non region of interest');
nroi_image=simg; 
logimg=imsubtract(orgimg,simg);
roi_image=logimg;
 subplot(2,1,2);
 imshow(logimg,[]), title('ROI-region of interest');

%%
%area calculation
[h,v]=size(im_bin);
tumour_area=(1/h)*(1/v)*(sum(sum(im_bin)));
display('area of tumour is :');
display (tumour_area);
D.tumour_area=tumour_area;
%%
%pixel calculation
[h,v]=size(im_bin);
number_pixel=(sum(sum(im_bin)));
display('number of pixel in the tumor area :');
display (number_pixel);
D.number_pixel=number_pixel;
%%
%subjective analysis
var1 = edge(im_bin,'canny');
var2 = imfuse(org_img,var1);
fuse_image=var2;
%figure, imshow(var2), title('subjective analysis');


if(p==1)
    ss=im2uint16(logimg);
else
    
    % for jpg images use im2uint8
    ss=im2uint8(logimg);
end


ss=im2uint8(ss);
seq=reshape(ss,[numel(ss) 1]);
img_hist=zeros(256,1);
for j=1:256
    img_hist(j)=sum(sum(ss==(j-1)));
end
count=img_hist';
count=count+1;
seq1=uint16(seq);
seq1=seq1+1;
code=arithenco(seq1,count);
decode=arithdeco(code,count,length(seq1));
decode=decode-1;
[r c]= size(ss);
img_out=reshape(decode,[r c]);
dc_roi=img_out;
%subplot(2,2,2);
 %imshow(img_out,[]),title('ROI after compression');
y=isequal(img_out,ss);
% display(y);
str3=sprintf('roi_before:roi_after = %f',y);
disp(str3);
%roi compression ratio
r_cr= (r*c*8)/(numel(code));
display(r_cr);
D.r_cr=r_cr;
%%
%roi bits per pixel
r_bpp=8/r_cr;
display(r_bpp);
D.r_bpp=r_bpp;
%%
if(p==1)
    
    % for dicom following conversion is needed
    s1=255/(max(max(simg)));
    simg=simg.*s1;
    simg=uint8(simg);

max_val = max(max(simg));
s1 = max_val/255;
simg = simg./s1;
simg=uint8(simg);

else
    simg=uint8(simg);
end

%%
imwrite(simg,'nroi_after.jpg','Mode','lossy','Quality',80);
n_info=imfinfo('nroi_after.jpg');
o=imread('nroi_after.jpg');
dc_nroi=o;
% subplot(2,2,4);
%imshow(o,[]),title('NROI after compression');
n_cr=(n_info.Height*n_info.Width)/n_info.FileSize;
% str1=sprintf('n_cr is %f',n_cr);
% disp(str1);
display(n_cr);
D.n_cr=n_cr;
n_bpp=8/n_cr;
% str2=sprintf('n_bpp is %f',n_bpp);
% disp(str2);
display(n_bpp);
D.n_bpp=n_bpp;
time=(r*c*8)/(28.4*1024);
D.Time=time;
str3=sprintf('Time taken is %f',time);
disp(str3);
display(time);
figure, imshow(var2), title('subjective analysis');
%%
%feature extraction


if(p==1)
    %     for DICOM
    logimg1= im2uint16(logimg);
    logimg1= uint8(logimg1);
    glcm1 = graycomatrix(logimg1);
    h_entropy = entropy(logimg1);
    
else
    % for jpeg
    glcm1 = graycomatrix(logimg);
    h_entropy = entropy(logimg);
end
%%
feature = graycoprops(glcm1,{'energy','contrast','homogeneity','correlation'});
display(feature);
D.h_entropy=h_entropy;
D.feature=feature;