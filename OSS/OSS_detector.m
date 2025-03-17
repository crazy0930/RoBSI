 function [ nkps, nM,M] = OSS_detector( im, sr, nk, K, s0)

F = 2 ^ ( 1 / 3 );
scale = 1 ./ F .^ ( 0:floor( log2( min( size( im ) ) ) ) - 2 - 1 );

[ H, W ] = size( im );
npts = round( H * W * K );%关键点数

Bw = sr + 4;

nkps = [  ];
nM = cell( length( scale ), 1 );
for i = 1:length( scale )
    if s0 == 0
        im2 = im;
    else
        sgi = s0 / scale( i );
        ga = fspecial( 'gaussian', [ 2 * round( 3 * sgi ) + 1, 2 * round( 3 * sgi ) + 1 ], sgi );
        im2 = imfilter( im, ga, 'replicate' );
    end
    if scale( i ) ~= 1
        im2 = imresize( im2, scale( i ) );
    end
    
    if size( im2, 1 ) < 30
        continue ;
    end
    
    [ U, M] = OSS_map( im2, sr, nk );
  
    %偏移子图构建
    U2 = U( Bw + 1:end  - Bw, Bw + 1:end  - Bw );
    [ r, c, rsp, csp ] = nonmaxsupptsgrid( U2, 4, 0.1, 1, 10000 );
%     
%     这行MATLAB代码调用了nonmaxsupptsgrid函数，该函数用于在图像中进行非
%     极大值抑制，以提取局部极大值点。该函数的输入参数包括：
% U2：一个二维矩阵，表示要在其上执行非极大值抑制的图像。
% 4：每个像素周围用于考虑非极大值抑制的局部邻域的大小。在这里，它设置为4。
% 0.1：用于确定像素是否为局部极大值的阈值值。亮度值小于此阈值的像素不被视为局部极大值并被抑制。这里的阈值设置为0.1。
% 1：图像中最大和最小强度值之间的比率。这用于归一化阈值值。
% 10000：要由函数返回的局部极大值点的最大数量。
% nonmaxsupptsgrid函数的输出包括：
% 
% r：一个向量，包含图像中局部极大值点的行坐标。
% c：一个向量，包含图像中局部极大值点的列坐标。
% rsp：每个局部极大值点的响应值（即强度值）。
% csp：每个局部极大值点的方向值。  
    
    if isempty( rsp )
        continue ;
    end
    kps = [ csp + Bw, rsp + Bw ];
    kps = kps ./ scale( i );
    
    nkps = [ nkps;[ kps, repmat( [ 1 ./ scale( i ), i ], size( kps, 1 ), 1 ) ], U2( sub2ind( size( U2 ), r, c ) ) ];
    
%     具体来说，代码中的 nkps 是一个矩阵，每一行代表一个关键点，包括该点的位置、尺度和朝向。
%     kps 是一个新的关键点集合，包含位置信息。scale 是一个数组，包含每个关键点的尺度。
%     U2 是一个与图像尺寸相同的二维数组，其中每个元素存储一个浮点数。r 和 c 是两个向量，
%     分别代表 kps 中每个关键点的行索引和列索引。

    
%     第一列是关键点的行索引，第二列是关键点的列索引
%     第一列是关键点的尺度信息，第二列是关键点的方向信息。
%     每一行代表一个关键点的响应值信息。
    nM{ i } = M; % M  48 特征描述子
end

if length( nkps( :, 5 ) ) > npts
    [ IV, ~ ] = sort( nkps( :, 5 ), 1, 'descend' );
    nkps = nkps( nkps( :, 5 ) >= IV( npts ), 1:4 );
    
% 这段 MATLAB 代码的作用是对矩阵 nkps 的第5列进行降序排列，
% 并将排序后的结果保存在变量 IV 中。另外，变量 ~ 表示一个占位符，用于忽略第二个返回值。
% 然后，这段代码将保留 nkps 矩阵中第5列的值大于等于 IV(npts) 的所有行，
% 并将结果保存在 nkps 矩阵的前4列中。其中 npts 是一个给定的整数索引值。
% 因此，这段代码可以用于对具有多列的矩阵进行排序并选择指定范围内的行。
else
    nkps = nkps( :, 1:4 );
end
end
