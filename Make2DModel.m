function pT = Make2DModel( P, TP, Mod )

x = P( :, 1 );
y = P( :, 2 );
Tp = TP.Tp;

switch Mod
    case 'Similarity'
        xT = ( Tp( 1 ) + Tp( 3 ) * x + Tp( 4 ) * y );
        yT = ( Tp( 2 ) - Tp( 4 ) * x + Tp( 3 ) * y );
        
    case 'Affine'
        xT = ( Tp( 1 ) + Tp( 2 ) * x + Tp( 3 ) * y );
        yT = ( Tp( 4 ) + Tp( 5 ) * x + Tp( 6 ) * y );
        
    case 'Bilinear'
        xT = ( Tp( 1 ) + Tp( 2 ) * x + Tp( 3 ) * y + Tp( 4 ) * x .* y );
        yT = ( Tp( 5 ) + Tp( 6 ) * x + Tp( 7 ) * y + Tp( 8 ) * x .* y );
        
    case 'Projective'
        xT = ( Tp( 1 ) + Tp( 2 ) * x + Tp( 3 ) * y ) ./ ( Tp( 7 ) * x + Tp( 8 ) * y + 1 );
        yT = ( Tp( 4 ) + Tp( 5 ) * x + Tp( 6 ) * y ) ./ ( Tp( 7 ) * x + Tp( 8 ) * y + 1 );
        
    case 'Polynomial'
        N = length( x );
        n = TP.n;
        ui = [ 3, 6, 10, 15, 21 ];
        u = ui( n );
        
        Polyb = [ ones( N, 1 ), x, y ...
            , x .* y, x .^ 2, y .^ 2 ...
            , ( x .^ 2 ) .* y, x .* ( y .^ 2 ), x .^ 3, y .^ 3 ...
            , ( x .^ 3 ) .* y, x .* ( y .^ 3 ), x .^ 4, y .^ 4, ( x .^ 2 ) .* ( y .^ 2 ) ...
            , ( x .^ 3 ) .* ( y .^ 2 ), ( x .^ 2 ) .* ( y .^ 3 ), x .^ 5, y .^ 5, ( x .^ 4 ) .* y, x .* ( y .^ 4 ) ];
        
        Ab = zeros( 2 * N, 42 );
        Ab( 1:2:end  - 1, 1:21 ) = Polyb;
        Ab( 2:2:end , 22:end  ) = Polyb;
        Ab( :, u + 1:21 ) = [  ];
        Ab( :, 2 * u + 1:end  ) = [  ];
        Lb = Ab * Tp;
        xT = ( Lb( 1:2:end  - 1 ) );
        yT = ( Lb( 2:2:end  ) );
        
    case 'LPiecewise'
        TRI = TP.TRI;
        T = pointLocation( TRI, [ x, y ] );
        xT = nan( size( x ) );yT = xT;
        
        for i = 1:length( x )
            if ~isnan( T( i ) )
                AfnPar = TP.Tp( T( i ), : );
                xT( i ) = AfnPar( 1 ) + AfnPar( 2 ) * x( i ) + AfnPar( 3 ) * y( i );
                yT( i ) = AfnPar( 4 ) + AfnPar( 5 ) * x( i ) + AfnPar( 6 ) * y( i );
            end
        end
        
    case 'MQuadrics'
        N = length( x );
        
        x1 = repmat( x, 1, length( TP.x ) );
        y1 = repmat( y, 1, length( TP.x ) );
        
        x2 = repmat( TP.x', N, 1 );
        y2 = repmat( TP.y', N, 1 );
        
        Ri = ( ( x1 - x2 ) .^ 2 + ( y1 - y2 ) .^ 2 + TP.d ) .^ .5;
        
        pT = Make2DModel( P, TP, 'Polynomial' );
        
        vx = Ri * TP.Tpx;
        vy = Ri * TP.Tpy;
        
        xT = pT( :, 1 ) + vx;
        yT = pT( :, 2 ) + vy;
        
    case 'TPS'
        N = length( x );
        
        x1 = repmat( x, 1, length( TP.x ) );
        y1 = repmat( y, 1, length( TP.x ) );
        
        x2 = repmat( TP.x', N, 1 );
        y2 = repmat( TP.y', N, 1 );
        
        Ri2 = ( x1 - x2 ) .^ 2 + ( y1 - y2 ) .^ 2 + TP.d;
        Ri = Ri2 .* log( Ri2 );
        Ri( isnan( Ri ) ) = 0;
        
        pT = Make2DModel( P, TP, 'Polynomial' );
        
        vx = Ri * TP.Tpx;
        vy = Ri * TP.Tpy;
        
        xT = pT( :, 1 ) + vx;
        yT = pT( :, 2 ) + vy;
        
    case 'MQuadricsI'
        N = length( x );
        
        x1 = repmat( x, 1, length( TP.x ) );
        y1 = repmat( y, 1, length( TP.x ) );
        
        x2 = repmat( TP.x', N, 1 );
        y2 = repmat( TP.y', N, 1 );
        
        Ri = exp(  - ( ( x1 - x2 ) .^ 2 + ( y1 - y2 ) .^ 2 + TP.d + eps ) );
        
        pT = Make2DModel( P, TP, 'Polynomial' );
        vx = Ri * TP.Tpx;
        vy = Ri * TP.Tpy;
        
        xT = pT( :, 1 ) + vx;
        yT = pT( :, 2 ) + vy;
    case 'TPSs'
        N = length( x );
        
        x1 = repmat( x, 1, length( TP.x ) );
        y1 = repmat( y, 1, length( TP.x ) );
        
        x2 = repmat( TP.x', N, 1 );
        y2 = repmat( TP.y', N, 1 );
        
        Ri2 = ( x1 - x2 ) .^ 2 + ( y1 - y2 ) .^ 2 + TP.d;
        Ri = Ri2 .* log( Ri2 );
        Ri( isnan( Ri ) ) = 0;
        
        pT = Make2DModel( P, TP, 'Polynomial' );
        vx = Ri * TP.Tpx;
        vy = Ri * TP.Tpy;
        
        xT = pT( :, 1 ) + vx;
        yT = pT( :, 2 ) + vy;
        
    case 'PW'
        pT = Make2DModel( P, TP, 'Polynomial' );
        N = length( x );
        
        x1 = repmat( pT( :, 1 ), 1, length( TP.x ) );
        y1 = repmat( pT( :, 2 ), 1, length( TP.x ) );
        
        x2 = repmat( TP.x', N, 1 );
        y2 = repmat( TP.y', N, 1 );
        Ri = ( ( x1 - x2 ) .^ 2 + ( y1 - y2 ) .^ 2 ) .^ .5;
        [ ~, Ind ] = sort( Ri, 2 );
        Ind = Ind( :, 1 );
        
        for i = 1:size( P, 1 )
            pT( i, : ) = Make2DModel( pT( i, : ), TP.Tpi{ Ind( i ) }, 'Polynomial' );
        end
        xT = pT( :, 1 );
        yT = pT( :, 2 );
        
    case 'PointWise'
        pT = Make2DModel( P, TP, 'Polynomial' );
        MdlKDT = KDTreeSearcher( [ TP.x, TP.y ] );
        [ KM, W ] = knnsearch( MdlKDT, pT, 'K', TP.k );
        
        dXi = TP.dX( KM );
        dYi = TP.dY( KM );
        W = ( 1 ./ ( W + eps ) .^ TP.pw );
        DX = sum( dXi .* W, 2 ) ./ sum( W, 2 );
        DY = sum( dYi .* W, 2 ) ./ sum( W, 2 );
        
        xT = pT( :, 1 ) - DX;
        yT = pT( :, 2 ) - DY;
        
    case 'WeightedMean'
        pT = Make2DModel( P, TP, 'Polynomial' );
        N = length( x );
        
        x1 = repmat( pT( :, 1 ), 1, length( TP.x ) );
        y1 = repmat( pT( :, 2 ), 1, length( TP.x ) );
        
        x2 = repmat( TP.x', N, 1 );
        y2 = repmat( TP.y', N, 1 );
        W = ( ( x1 - x2 ) .^ 2 + ( y1 - y2 ) .^ 2 ) .^ .5 + eps;
        W = 1 ./ W .^ TP.pw;
        
        DX = zeros( N, 1 );DY = DX;
        for i = 1:N
            DX( i ) = sum( TP.dX .* W( i, : )' ) / sum( W( i, : ) );
            DY( i ) = sum( TP.dY .* W( i, : )' ) / sum( W( i, : ) );
        end
        xT = pT( :, 1 ) - DX;
        yT = pT( :, 2 ) - DY;
        
end

pT = [ xT, yT ];

