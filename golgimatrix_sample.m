
function golgimatrixwin12
%Copyright, Dmitry Yampolsky 2022
%golgimatrix project

%code sample

%Kai M. Bracey, Kung-Hsien Ho, Dmitry Yampolsky, Guogiang Gu, Irina Kaverina, William R. Holmes,
%Microtubules Regulate Localization and Availability of Insulin Granules in Pancreatic Beta Cells,
%Biophysical Journal, Volume 118, Issue 1, 2020, Pages 193-206

%dependencies:

% CircStat2012a
% SPHERE_VORONOI
% resize
% GridSphere
% UniformSampling
% imgaussian

clear all;

noisedSwitch = false;
generateMswitch = true;
ncases = 1;
matrixsize = 2^2;
nmatrixes = 5;

drawSwitch = false;

npixels = 1000;
binCeil = 314;

  path(path,'./CircStat2012a');
  path(path,'./SPHERE_VORONOI');
  path(path, './resize');
  path(path,'./GridSphere');
  path(path,'./UniformSampling');
  

    %load data
    [filename, pathname]=  uigetfile('*.gif');
               dataStruct = importdata([pathname,filename]);
               cdata=dataStruct.cdata;
                  cdata = squeeze(cdata);
                  
                  if ndims(cdata)~=3
                       error('Image file must be a stack');
                  end
                  
                  if length(unique(cdata))~=2 %
    
                      cdata=double(cdata);
                       cdata2 = resize(double(cdata),[100,100,size(cdata,3)]);
                      
                      threshseek=2;
                      threshmark = -1;
               
                     pixdiff =-3;
                     pixdiff2=-3;
                     
                     while(abs(pixdiff)>2)
                           while sign(pixdiff2)==sign(pixdiff)
                               pixdiff=pixdiff2;
                              threshseek = threshseek +threshmark;
                               cdata2 =  irish_threshold(cdata,100,threshseek);
                               pixdiff2 = (sum(cdata2(:))-npixels);
                           end
                           
                         threshmark = -threshmark*.5;
                           pixdiff=pixdiff2;
                     end
                     
                     cdata=cdata2;
                     
                  else
                  cdata(:) = cdata(:)/max(cdata(:));%make binary
                  
                  if  numel(cdata)/2 < sum(cdata(:))%if background(majority) is black
                                                           %invert
                     cdata(:)=~cdata(:);                  
                  end
                  end
        %image data is now binary
        
        if size(cdata,1)~=size(cdata,2)%make image size square
                                                   %(image is presumably
                                                   %approximstely square)
              cdata = make_square(cdata);
        end
      cdataOriginal = cdata;
      %end load data   
      
                  %trans to vector
                threeDmatrixcoords=  make_spherical_vectors(cdata);
                  
                   resultsarray = [];   
                   

     uniformmodel = ParticleSampleSphere('N',40);
        [IOD,Dhist] =  IOD_from_norm(threeDmatrixcoords,uniformmodel);
        
        
       maxfreq = 314;
        maxfreq = 50;
      Dhistceiled=min(Dhist,maxfreq); 
        
        disp('IOD = ');
        disp(var(Dhistceiled)/mean(Dhistceiled));
        figure;
         plot3(threeDmatrixcoords(:,1),threeDmatrixcoords(:,2),threeDmatrixcoords(:,3),'b.');
        return
        
        
    function pl3(in)
        plot3(in(:,1),in(:,2),in(:,3),'o');
          
        
    function  vecout = make_spherical_vectors(binaryin)
        
         [x,y,z]=meshgrid(0:size(binaryin,1)-1,0:size(binaryin,2)-1,0:size(binaryin,3)-1);
                  
                 % 0:1 -> -1:1
                  x=2*(x/size(binaryin,1))-1;
                  y=2*(y/size(binaryin,2))-1;
                  z=2*(z/size(binaryin,3))-1;
                 % 
                  %make 3d coord array
                vectmp= [x(logical(binaryin)),y(logical(binaryin)),z(logical(binaryin))];
              %normalize to unit sphere
           vecout = vectmp ./repmat(sum(vectmp.^2,2).^.5,[1,3]);
                    
        
    function dataout = make_square(datain)%transform image's first 2 dimensions to be equal size

                             [~,dimTmp] = max([size(datain,1),size(datain,2)]);
                      diffTmp = abs(size(datain,1)-size(datain,2));
                      diffTmp1 = ceil(diffTmp/2);
                      diffTmp2 = floor(diffTmp/2);
                      
                      cdata2 = zeros([size(datain,dimTmp),size(datain,dimTmp),size(datain,3)]);
                      if dimTmp == 2
                          cdata2((diffTmp1+1):end-diffTmp2,:,:)=datain;
                      else%dimTmp == 2
                          cdata2(:,(diffTmp1+1):end-diffTmp2,:)=datain;
                      end
                      dataout = cdata2;
        


function [xout, yout, zout]=geodesic_grid(scalein,randfactorin)
                                                    %in degrees
                                                    
[latGridInDegrees, longGridInDegrees] = GridSphere(scalein);
[xout yout zout]=sph2cart(degtorad(longGridInDegrees),degtorad(latGridInDegrees),ones(length(longGridInDegrees),1));


    function vecout = make_gauss_model(n_in,res_in,sigma_in,r_in,thresh_in)        
 uniformmodel_gauss = ParticleSampleSphere('N',n_in);
   uniformmodel_gauss=(uniformmodel_gauss+1)/2;
   uniformmodel_gauss=floor(uniformmodel_gauss*res_in)+1;%cube
 zeromesh=zeros(res_in,res_in,res_in);
 
 for i = 1:n_in
 zeromesh(uniformmodel_gauss(i,1),uniformmodel_gauss(i,2),uniformmodel_gauss(i,3))=1;
 end
 
 uniformmodel_gauss2=imgaussian(zeromesh,sigma_in,r_in);
 uniformmodel_gauss2_bin = (uniformmodel_gauss2>thresh_in);
  vecout = make_spherical_vectors(uniformmodel_gauss2_bin);
  geo=geodesic_grid(length(vecout),0);
 
function ptsout = random_distribution_model(npts)

    %a random data model  
  ex=2*rand(npts,3)-1;%-1:1 cube of random points
  ex2=ex(((ex(:,1).^2+ex(:,2).^2+ex(:,3).^2)<1),:);%cutout the unit sphere
  ptsout=ex2./repmat(sqrt(sum(ex2.^2,2)),[1,3]);%move the sphere to the sphere surface by normalizing

