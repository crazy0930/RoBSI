function  [Multi_gradient,cmpc,dMSOG] = Dense_MSFG_Features(image,div)
%%  采用四方向扩展 Sobel算子模板 [0, 45, 90, 135]
%% Multi-directed two-order gradients based on 3D Gaussian filtering
Gimage=im2double(image);
[M,N]=size(image);
P=div;
%% 初始化
Multi_gradient=zeros(M,N,P);
Multi_angle=zeros(M,N,P);

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
    or1   = atan2( gradient_y_1, gradient_x_1);
    Multi_angle(:,:,i)=or1;
    Multi_gradient(:,:,i)=gradient_1;
%     figure, imshow(Multi_gradient(:,:,j),[]);
%     figure, imshow(Multi_angle(:,:,j),[]);
end

for j=1:1:P-4
    gradient_x_1=imfilter(Gimage,H{j},'replicate');
    gradient_y_1=imfilter(Gimage,H{j+1},'replicate');
    gradient_1=sqrt(gradient_x_1.^2+gradient_y_1.^2);
    %First order Sobel filter
    gradient_x_2=imfilter(gradient_1,H{j},'replicate');
    gradient_y_2=imfilter(gradient_1,H{j+1},'replicate');
    gradient_2=sqrt(gradient_x_2.^2+gradient_y_2.^2);
    %Second order Sobel filter
    or   = atan2( gradient_y_2, gradient_x_2); 
    Multi_angle(:,:,j+4)=or;
    Multi_gradient(:,:,j+4)=gradient_2;
%     figure, imshow(Multi_gradient(:,:,j),[]);
%     figure, imshow(Multi_angle(:,:,j),[]);
end
cmpc= normalize(Multi_angle, P);
% figure, imshow(cmpc(:,:,j),[]);
sigma =0.8;
f = fspecial( 'gaussian', max( 1, fix( 6 * sigma + 1 ) ), sigma );
f = cat( 3, 1 * f, 3 * f, 1 * f );
dMSOG1 = convn(Multi_gradient, f, 'same' );
sum1 = sum( dMSOG1, 3 );
dMSOG = dMSOG1 ./ ( sum1 + 0.000000001 );
% figure, imshow(dMSOG(:,:,j),[]);
end

%%  进行多通道三维卷积
function  cmpc= normalize(temp, div)
g = fspecial('gaussian',3,0.5);%Two-dimensional Gaussian kernel
o = div;
for i = 1:o
    temp(:,:,i) = conv2(temp(:,:,i), g, 'same'); 
end

[r,c] = size(temp(:,:,1));
% increase all zero layers
temp_array = zeros(r,c,o+2);
for i = 2:o+1
    temp_array(:, :, i) = temp(:, :, i - 1);
end

delta = 0.0000000001;
denominator_normalized = zeros(r,c,1) + delta;

% z-axis filtering by a kernel = [1,3,1]
cmpc = zeros(r,c,o);
for i = 1:o
    cmpc(:, :, i) = temp_array(:, :, i) + 3 * temp_array(:, :, i + 1) + temp_array(:, :, i + 2);
    denominator_normalized = denominator_normalized + cmpc(:, :, i).^(2);
end

% L2 norm regularization
denominator_normalized = denominator_normalized.^(0.5);
for i = 1:o
    cmpc(:, : , i) = cmpc(:, : , i)./denominator_normalized;
end
end


