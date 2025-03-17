function [ U, imN] = OSS_map( im, sr, nk )%U自相似特征图
im = im * 255;
%构建子图
[ R, C ] = size( im );
RP = R - 2 * sr;
CP = C - 2 * sr;
SP = RP * CP;

sr_mask = strel( 'disk', sr, 0 ).Neighborhood;
sr_mask( sr + 1, sr + 1 ) = false;
sr_idx = find( sr_mask );

L = length( sr_idx );

h = fspecial( 'disk', 2 );%因此，代码h = fspecial( 'disk', 2 );将创建一个半径为2的圆形滤波器，并将其存储在名为'h'的变量中，以便进一步使用。这个圆形滤波器可以用于图像处理中的各种滤波任务，例如平滑、锐化或边缘检测等。

imC = im( 1 + sr:end  - sr, 1 + sr:end  - sr );
imN = single( zeros( RP, CP, L ) );

[ pn, qn ] = ind2sub( [ size( sr_mask, 1 ), size( sr_mask, 2 ) ], sr_idx );
% 具体来说，[size(sr_mask,1), size(sr_mask,2)]返回sr_mask的行数和列数的大小，
% 创建一个2D数组作为ind2sub函数的第一个输入参数。sr_idx是一个线性索引，
% 表示sr_mask中某个元素在一个按列堆叠的向量中的索引位置。
% ind2sub函数的作用是将这个线性索引转换为矩阵中的行列索引，
% pn和qn分别对应于sr_idx在矩阵中的行号和列号。
% 因此，这行代码的结果是将sr_idx转换为pn和qn，这些变量可以用于定位sr_mask矩阵中的元素。

for i = 1:L / 2
    p = pn( i );q = qn( i );
    imT = im( p:p + RP - 1, q:q + CP - 1 );
    imT = abs( imT - imC );
    imT = imfilter( imT, h, 'replicate' );
    imN( :, :, i ) = circshift( imT, [ p - 1 - sr, q - 1 - sr ] );
    imN( :, :, L - i + 1 ) = imT;
end

% 
% 这段MATLAB代码中，L是一个偶数，pn和qn是先前计算出的行列索引数组。这段代码的目的是计算每个参考块周围的非参考块。
% 具体来说，该循环迭代了L/2次（即块数的一半），在每次迭代中，首先获取参考块的行列索引（p和q）。
% 接下来，从原始图像中提取大小为RP×CP的图像块，并将其减去参考块。
% 然后，通过使用卷积核h对差异图像进行平滑处理，即imT = imfilter( imT, h, 'replicate' )。
% 接下来，使用circshift函数将非参考块平移并存储在imN中。其中第一个存储的非参考块是距离
% 参考块最近的非参考块（相对于参考块的左上角），而第二个存储的非参考块则是距离参考块最远的非参考块（相对于参考块的右下角）。
% 最终，该循环会计算出与每个参考块对应的非参考块，并将这些非参考块存储在imN中，imN的大小为[RP×CP×L/2]。

clear imC imT;

imN = reshape( imN, SP, L );
%使用reshape函数将imN重新塑造为一个大小为[SP×L]的矩阵，其中每一列代表一个非参考块。

[ IV, ~ ] = sort( imN, 2 );
%48个小块
% 最终，sort函数将返回一个大小与imN相同的矩阵IV，其中每一列包含按升序
% 排列的非参考块元素。例如，IV的第一列包含最相似的非参考块，而IV的最后一列包含最不相似的非参考块。
U = sum( IV( :, 1:nk ), 2 ) ./ nk;

% 这行MATLAB代码的作用是计算每个像素的非参考块的加权平均值，并将结果存储在U中。
% 具体来说，IV是一个大小为[SP×L]的矩阵，其中每一列包含已排序的非参考块元素。
% nk是一个整数，表示我们使用IV中的前nk个元素来计算加权平均值。因此，IV( :, 1:nk )
% 是一个大小为[SP×nk]的矩阵，包含IV中前nk个非参考块元素。对这个矩阵的每一行求和，将返回一个
% 大小为[SP×1]的列向量，表示每个像素的非参考块加权平均值。
% 最后，我们将这些加权平均值除以nk，以获得每个像素的非参考块的平均值。这些值将被存储在列向量U中。

U = reshape( U, RP, CP );
U = padarray( U, [ sr, sr ], 'replicate', 'both' );
% 代码的含义是将名为"U"的数组在顶部，底部，左侧和右侧各添加了 sr 行和 sr 列，
% 并使用'replicate'填充方法来填充新添加的元素。'both'参数指定在所有四个方向上进行填充。
imN = reshape( imN, RP, CP, L );
imN = padarray( imN, [ sr, sr, 0 ], 'replicate', 'both' );

end
