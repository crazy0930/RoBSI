function  [Multi_gradient,dMSOG] = Dense_MFOG_Features2(image,div)
%%  采用四方向扩展 Sobel算子模板 [0, 45, 90, 135]
%% Multi-directed two-order gradients based on 3D Gaussian filtering
Gimage=im2double(image);
[M,N]=size(image);
P=div;
%% 初始化
Multi_gradient=zeros(M,N,P);
h1=[1,2,1;0,0,0;-1,-2,-1];            % Sobel  0
h2=[2,1,0;1,0,-1;0,-1,-2];            % Sobel  45
h3=[1,0,-1;2,0,-2;1,0,-1];            % Sobel  90
h4=[0,-1,-2;1,0,-1;2,1,0];            % Sobel  135
h5=h1;
H={h1,h2,h3,h4,h5};

%%  计算多向Sobel算子结果
for i=1:1:P-4
    %First order Sobel filter
    gradient_x_1=imfilter(Gimage,H{i},'replicate');
    gradient_y_1=imfilter(Gimage,H{i+1},'replicate');
    gradient_1=sqrt(gradient_x_1.^2+gradient_y_1.^2);
    Multi_gradient(:,:,i)=gradient_1;
%     figure, imshow(Multi_gradient(:,:,i),[]);
end
for j=1:1:P-4
    gradient_x_1=imfilter(Gimage,(H{j}),'symmetric','conv');
    gradient_y_1=imfilter(Gimage,(H{j+1}),'symmetric','conv');
    gradient_2=sqrt(gradient_x_1.^2+gradient_y_1.^2);
    Multi_gradient(:,:,j+4)=gradient_2;
%     Ng=mapminmax(gradient_2,0,1);
%     map=addcolorplus(335);
%     figure, imshow(Ng, 'Colormap' ,map);
end

sigma =0.7;
f = fspecial( 'gaussian', max( 1, fix( 6 * sigma + 1 ) ), sigma );
f = cat( 3, 1 * f, 3 * f, 1 * f );
dMSOG1 = convn(Multi_gradient, f, 'same' );
sum1 = sum( dMSOG1, 3 );
dMSOG = dMSOG1 ./ ( sum1 + 0.000000001 );
%  map=addcolorplus(339);
%  figure, imshow(dMSOG(:,:,j), 'Colormap' ,map);
end
