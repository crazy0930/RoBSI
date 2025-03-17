function des = OSS_desG( M, KPS, PS, o, br, ba, sr, sr_mask )
n = ( br - 1 ) * ba + 1;
des = zeros( n * o, size( KPS, 2 ) );
ign = zeros( 1, size( KPS, 2 ) );

for k = 1:size( KPS, 2 )
    x = round( KPS( 1, k ) );
    y = round( KPS( 2, k ) );
    a = KPS( 4, k );
    hs = round( PS / 2 );
    
    x1 = max( 1, x - ( hs ) );
    y1 = max( 1, y - ( hs ) );
    x2 = min( x + ( hs ), size( M, 2 ) );
    y2 = min( y + ( hs ), size( M, 1 ) );
    
    if y2 - y1 ~= 2 * hs || x2 - x1 ~= 2 * hs
        ign( k ) = 1;continue ;
    end
    
    [ log_angle1, log_angle2, log_amp2 ] = OSS_mask( sr, hs, 2 * o, br, ba, a );
    M2 = M( y1:y2, x1:x2, : );
    M8 = zeros( PS + 1, PS + 1, o );
    for i = 1:o
        idx = sr_mask( log_angle1 == i );
        M8( :, :, i ) = mean( M2( :, :, idx ), 3 );
    end
    [ ~, patch ] = min( M8, [  ], 3 );
    
    
    tmp_des = zeros( n, o );
    mm = 1;
    idx = find( log_amp2 == 1 );
    clip = patch( idx );
    hv = hist( clip( clip > 0 ), 1:o );
    tmp_des( mm, : ) = hv( 1:o ) ./ length( idx );
    mm = mm + 1;
    for i = 2:br
        for j = 1:ba
            idx = find( log_amp2 == i & log_angle2 == j );
            clip = patch( idx );
            hv = hist( clip( clip > 0 ), 1:o );
            tmp_des( mm, : ) = hv( 1:o ) ./ length( idx );
            mm = mm + 1;
        end
    end
    tmp_des = tmp_des( : );
    
    if norm( tmp_des ) ~= 0
        tmp_des = tmp_des / norm( tmp_des );
    end
    des( :, k ) = tmp_des;
end
des = struct( 'kps', KPS( :, ign == 0 )', 'des', des( :, ign == 0 )' );
end
