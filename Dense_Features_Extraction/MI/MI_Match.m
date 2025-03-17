function  [CP_Ref ,CP_Sen] =MI_Match(im_Ref,im_Sen,Detectors,KPT_NUMS,H)
tic
im_Ref = im2double(im_Ref);
im_Sen = im2double(im_Sen);
[~,~,k3] = size(im_Ref);
if k3 > 1
    im_Ref = rgb2gray(double(im_Ref));
end
im_Ref = double(im_Ref);
[~,~,k4] = size(im_Sen);
if k4 > 1
    im_Sen = rgb2gray(double(im_Sen));
end
im_Sen = double(im_Sen);

[im1height,im1width] = size(im_Ref);
%% 分块提取特征点
corrRad = 50;                     %模版半径
rad = 10;                            %搜索半径，默认：10
marg=corrRad+rad+2;      %边界
%提取分块的harris特征点
im1 = im_Ref(marg:im1height-marg,marg:im1width-marg);
% Value = harrisValue(im1);
% [r,c,~,~] = nonmaxsupptsgrid(Value,3,0.001,10,3);
% points1 =[r,c] + marg - 1;
% m1_points = detectFASTFeatures(im11,'MinContrast',0.05);
% m1_points = m1_points.selectStrongest(500); 
% points1 =[round(m1_points.Location(:,1)),round(m1_points.Location(:,2))] + marg - 1;
switch(Detectors)
    case 'B_Harris'   %  多向双层一阶梯度卷积
        block_kpts = 10;
        B_nums = KPT_NUMS/(block_kpts*block_kpts);
        im1 = im1*255;
        Value = harrisValue(im1);
        [r,c,rsubp,cubp] = nonmaxsupptsgrid(Value,3,0.3,block_kpts,B_nums); 
        points1 =[r,c] + marg - 1;
     case 'ANMS'   %  多向双层一阶梯度卷积
         a=max(im1(:));  b=min(im1(:));  imn=(im1-b)/(a-b); 
         m1_points= ANMS(imn,KPT_NUMS);
         points1 =[round(m1_points.Location(:,2)),round(m1_points.Location(:,1))] + marg - 1;
%%  提出的新特征点提取
    case 'GBMS'   %  多向双层一阶梯度卷积
        kpts=GridsBoxPoints(im1,KPT_NUMS,0.2);
        r =kpts(:,1);   c =kpts(:,2);
        points1 =[r,c] + marg - 1;
    case 'FAST'   %  多向双层一阶梯度卷积
        m1_points = detectFASTFeatures(im1,'MinContrast',0.05);
        m1_points = m1_points.selectStrongest(KPT_NUMS); 
        points1 =[round(m1_points.Location(:,2)),round(m1_points.Location(:,1))] + marg - 1;
end

%% 计算影像梯度特征
% [g1,or1,~] =imgrad(im_Ref);
% [g2,or2,~] =imgrad(im_Sen);

pNum = size(points1,1);
disthre = 5; %error threshold   1.5
C = 0;%the number of correct match 
CM = 0 ;%the number of total match 
C_e = 0;%the number of mismatch

for n = 1: pNum
    
     %the x and y coordinates in the reference image
     X_Ref=points1(n,2);
     Y_Ref=points1(n,1);
   
    inputData1 = (im_Ref(Y_Ref-corrRad:Y_Ref+corrRad,X_Ref-corrRad:X_Ref+corrRad));
    %transform the (x,y) of reference image to sensed image by the geometric relationship of check points 
%     tempCo22 = [X_Ref,Y_Ref];

    tempCo = [X_Ref;Y_Ref;1];
    tempCo1 = H*tempCo;
    
    %tranformed coordinate (X_Sen_c, Y_Sen_c)
    X_Sen_c = tempCo1(1);
    Y_Sen_c = tempCo1(2);
    X_Sen_c1=round(tempCo1(1));
    Y_Sen_c1 =round(tempCo1(2)); 
    %judge whether the transformed GCP is out the boundary of right image.

    if (X_Sen_c1 < marg+1 | X_Sen_c1 > size(im_Sen,2)-marg | Y_Sen_c1<marg+1 | Y_Sen_c1 > size(im_Sen,1)-marg)
        %if out the boundary, this produre enter the next cycle
        continue;
    end
    corr = zeros(2*rad + 1);
    for i = -rad:rad
        for j = -rad:rad
            Y_Sen_c2 = Y_Sen_c1+i;
            X_Sen_c2 = X_Sen_c1+j;
            inputData2 = (im_Sen( Y_Sen_c2-corrRad: Y_Sen_c2+corrRad,X_Sen_c2-corrRad:X_Sen_c2+corrRad));
            temp=MI(inputData1,inputData2,127); %calculate MI
            corr(i + rad +1,j + rad + 1)=temp;
%             figure,imshow(corr,[]);
        end
    end
     maxCorr = max(max(corr));
     max_index = find(corr == maxCorr);
      if(size(max_index,1) > 1)
            % if two maxumal appear, it go to nexe cycle;
         continue;
      end
      [max_i,max_j] = ind2sub(size(corr),max_index);
       %the (matchY,matchX) coordinates of match
      matchY = Y_Sen_c1-rad + max_i -1;
      matchX = X_Sen_c1-rad + max_j -1;
      % calculate the match errors         
       diffY = abs(matchY- Y_Sen_c );
       diffX = abs(matchX- X_Sen_c );

       diff = sqrt(diffX.^2+diffY.^2);
      % calculate the numbers of correct match, mismatch and total match
      if diff <= disthre
          C = C+1;
          corrp(C,:)=[X_Ref,Y_Ref,matchX,matchY,diff];
      else
          C_e = C_e+1;
          corrp_e(C_e,:)=[X_Ref,Y_Ref,matchX,matchY,diff];
      end
      CM = CM + 1;
end
%the correct match ratio
precision = C/CM;
% x= sprintf('the correct match ratio (i.e., precision) is %4.3f',precision);
% disp(x)
fprintf('the total matching time is %f\n',toc);
CP_Ref = corrp(:,1:2);
CP_Sen = corrp(:,3:4);
end

function [g1,or1,or2] =imgrad(im1)
    h1=[1,2,1;0,0,0;-1,-2,-1];            % Sobel  0
    gradient_x_1=imfilter(im1, h1,'replicate');
    gradient_y_1=imfilter(im1, h1','replicate');
    g1 =sqrt(gradient_x_1.^2+gradient_y_1.^2);
    or1 = atan2(- gradient_y_1, gradient_x_1 );
    or2 = atan( gradient_y_1 ./ ( gradient_x_1 + 0.00000001 ) );
end