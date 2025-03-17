function dCSOG = Dense_CSOG( im, orbin, sigma)
if nargin < 3
    orbin = 9;
    sigma = 0.8;
end
[ h, w, o ] = size( im );
if o == 3
    im = rgb2gray( im );
end
orhistbig = zeros( h, w, orbin + 2 );
[g,~,or1] =imgrad(im);
or1( or1 < 0 ) = or1( or1 < 0 ) + pi;
theta = pi / orbin;
ortemp = ( or1 + theta / 2 ) / theta + 1;
orInt = floor( ortemp );
orInt1 = orInt + 1;
orFrac = ortemp - orInt;
orInt_val = g .* ( 1 - orFrac );
orInt1_val = g .* orFrac;
for i = 1:h
    for j = 1:w
        orhistbig( i, j, orInt( i, j ) ) = orInt_val( i, j );
        orhistbig( i, j, orInt1( i, j ) ) = orInt1_val( i, j );
    end
end
f = fspecial( 'gaussian', max( 1, fix( 6 * sigma + 1 ) ), sigma );
f = cat( 3, 1 * f, 3 * f, 1 * f );
dCSOG1 = convn( orhistbig, f, 'same' );
dCSOG = dCSOG1( :, :, 2:end  - 1 );
sum1 = sum( dCSOG, 3 );
dCSOG = dCSOG ./ ( sum1 + 0.000000001 );
end

function [g,or,or1] =imgrad(im1)
    h1=[1,2,1;0,0,0;-1,-2,-1];            % Sobel  0
    gradient_x_1=imfilter(im1, h1,'replicate');
    gradient_y_1=imfilter(im1, h1','replicate');
    g1 =sqrt(gradient_x_1.^2+gradient_y_1.^2);
    gradient_x_2=imfilter(g1, h1,'replicate');
    gradient_y_2=imfilter(g1, h1','replicate');
    g =sqrt(gradient_x_2.^2+gradient_y_2.^2);
    or = atan2(- gradient_y_2, gradient_x_2 );
    or1 = atan( gradient_y_2 ./ ( gradient_x_2 + 0.00000001 ) );
end
