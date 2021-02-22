close all
clear all
maxproj_dir='D:\2020-12-01-smFISH\StandardCulture_FISH_dapiNodalLefty\MIP_dapi_nodal_lefty\';
ff = dir(maxproj_dir);
segm_chan = 1;% dapi
q_chan =2;% [dapi nodal lefty] [1 2 3]
chan_nm = {'Dapi','Nodal','Lefty'};
cd(maxproj_dir);
mask=[];
fnm = struct;
ilastik = struct;
prob_thresh=0.9;
% params for kymograph for nodal movie
pxl_to_micron = 0.3250;%( pxltomicron = 0.325 - 40X)
small_stuff=500;
mean_expr=[];
v = [0.5 1 2 4 8];%
v = 1;
data = struct;
mask = struct;
save_nm = [maxproj_dir chan_nm{q_chan} '_Mtsr_data.mat'];
for tp=1:size(v,2)   
q = 1;
for jj=1:size(ff,1)
    if v(tp) >=1
    test_str = ['mTeSR'];%'Activin' num2str(v(tp)) 'h_MIP'
    else
    test_str = ['mTeSR'];  %  'Activin30m_MIP'
    end
    
    if  ~isdir(ff(jj).name) &&  ~isempty(regexp(ff(jj).name,[test_str ],'ONCE')) && ~isempty(regexp(ff(jj).name,['_w000' num2str(segm_chan-1) '.tif'],'ONCE')) %'mTeSR' 
      
        fnm(q).name =  ff(jj).name;        
        ilastik(q).name = [fnm(q).name(1:end-4) '_Probabilities.h5'];% dapi
        %disp(ilastik(q).name);
        mask(q).all = imfill(bwareaopen(readIlastikProbMask(ilastik(q).name,prob_thresh),small_stuff),'holes');                 
        disp([ff(jj).name(1:size(ff(jj).name,2)-5) num2str(q_chan-1) '.tif']);
        img = imread([ff(jj).name(1:size(ff(jj).name,2)-5) num2str(q_chan-1) '.tif']);% stain
        img_no_bg = img-mean(img(~imdilate(mask(q).all,strel('disk',30))));
        imshow(img_no_bg,[]);
        %imshowpair(img,mask(q).all);
        stats = regionprops(imdilate(mask(q).all,strel('disk',10)),img_no_bg,'Area','Centroid','MeanIntensity');
        mean_expr(q,1) = mean(cat(1,stats.MeanIntensity));
        data(tp).expr = mean_expr;       
        q = q+1;
    end
   
end
end
close all
for ii=1:size(v,2)
figure(1), errorbar(v(ii),mean(cat(1,data(ii).expr)),std(cat(1,data(ii).expr)),'-rp','LineWidth',2);hold on
end
xlabel('Time, hrs')
ylabel('Mean pixel intensity in dilated nuclei')
title(chan_nm{q_chan})
h6 = figure(1);
h6.CurrentAxes.FontSize = 12;
h6.CurrentAxes.LineWidth = 2;
xlim([0 10])

%save(save_nm,'maxproj_dir','segm_chan','q_chan','v','data','chan_nm','small_stuff');

