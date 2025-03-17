function kps2 = OSS_orient( M, kps, sr, nb, sr_mask )
T = 0.8;
if sr == 4
    n_grid = [ 4, 0;3, 1;3, 2;2, 3;1, 3;
        0, 4; - 1, 3; - 2, 3; - 3, 2; - 3, 1;
        - 4, 0; - 3,  - 1; - 3,  - 2; - 2,  - 3; - 1,  - 3;
        0,  - 4;1,  - 3;2,  - 3;3,  - 2;3,  - 1;4, 0 ];
elseif sr == 5
    n_grid = [ 5, 0;4, 1;4, 2;4, 3;3, 3;3, 4;2, 4;1, 4;
        0, 5; - 1, 4; - 2, 4; - 3, 4; - 3, 3; - 4, 3; - 4, 2; - 4, 1;
        - 5, 0; - 4,  - 1; - 4,  - 2; - 4,  - 3; - 3,  - 3; - 3,  - 4; - 2,  - 4; - 1,  - 4;
        0,  - 5;1,  - 4;2,  - 4;3,  - 4;3,  - 3;4,  - 3;4,  - 2;4,  - 1;5, 0 ];
else
    n_grid = [ 3, 0;2, 1;2, 2;1, 2;
        0, 3; - 1, 2; - 2, 2; - 2, 1;
        - 3, 0; - 2,  - 1; - 2,  - 2; - 1,  - 2;
        0,  - 3;1,  - 2;2,  - 2;2,  - 1;3, 0 ];
end
%边缘点角度
n_angle = atan2( n_grid( :, 2 ), n_grid( :, 1 ) ) * 180 / pi;

n_angle( n_angle < 0 ) = n_angle( n_angle < 0 ) + 360;
n_angle( end  ) = 360;

kpn = 0;
kps2 = zeros( size( kps, 1 ), 4 );

n = nb;
for i = 1:size( kps, 1 )%特征点集合
    c0 = round( kps( i, 1 ) );%图像中的位置
    r0 = round( kps( i, 2 ) );
    %根据特征点位置计算出一个局部区域，并将其转换成线性索引。
    idx = sr_mask( sub2ind( size( sr_mask ), n_grid( :, 2 ) + sr + 1, n_grid( :, 1 ) + sr + 1 ) );
    idx2 = sub2ind( size( M ), repmat( r0, length( idx ), 1 ), repmat( c0, length( idx ), 1 ), idx );
    n_bins = M( idx2 );%第6行根据线性索引，从一个预先计算好的直方图 M 中获取一组直方图值。
    bins = interp1( n_angle, n_bins, 0:360 / nb:360 - 360 / nb, 'pchip' );%将这组直方图值插值成一个 36个 bin 的直方图。
    hist = ( max( bins ) - bins ) ./ ( max( bins ) - min( bins ) );
    %计算出一个归一化的直方图，用于描述当前特征点的局部图像特征
    [ val, idx ] = sort( hist, 'descend' );
    val_th = val( 1 ) * T;
    count = 0;
    for j = 1:length( idx )
        kk = idx( j );
        if kk == 1
            k1 = nb;k2 = 2;
        elseif kk == nb
            k1 = nb - 1;k2 = 1;
        else
            k1 = kk - 1;k2 = kk + 1;
        end
        if ( hist( kk ) > hist( k1 ) && hist( kk ) > hist( k2 ) && hist( kk ) > val_th )
            bin = kk - 1 + 0.5 * ( hist( k1 ) - hist( k2 ) ) / ( hist( k1 ) + hist( k2 ) - 2 * hist( kk ) );
            if ( bin < 0 )
                bin = nb + bin;
            elseif ( bin >= nb )
                bin = bin - nb;
            end
            kpn = kpn + 1;
            kps2( kpn, 1 ) = kps( i, 1 );
            kps2( kpn, 2 ) = kps( i, 2 );
            kps2( kpn, 3 ) = kps( i, 3 );
            kps2( kpn, 4 ) = ( 360 / nb ) * bin;
            count = count + 1;
        end
        if count == 2
            break ;
        end
    end
end
end
% 对这个归一化的直方图进行一些处理，具体包括对直方图的排序、
% 计算阈值等等。如果符合一定条件，则将当前特征点及其对应的直方
% 图 bin 保存到 kps2 变量中，这样就得到了一个特征点的描述向量。