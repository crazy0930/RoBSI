function [ log_angle1, log_angle2, log_amp2 ] = OSS_mask( R1, R2, d1, br, ba, a )
cos_t = cos(  - a / 180 * pi + pi / 2 );
sin_t = sin(  - a / 180 * pi + pi / 2 );

d = d1;
X =  - R1:1:R1;
Y =  - R1:1:R1;
[ XX, YY ] = meshgrid( X, Y );
c_rot = XX * cos_t - YY * sin_t;
r_rot = XX * sin_t + YY * cos_t;

log_angle = atan2( r_rot, c_rot );
log_angle = log_angle / pi * 180;
log_angle( log_angle < 0 ) = log_angle( log_angle < 0 ) + 360;
log_amp = log2( sqrt( c_rot .^ 2 + r_rot .^ 2 ) );

log_angle = round( log_angle * d / 360 );
log_angle( log_angle <= 0 ) = log_angle( log_angle <= 0 ) + d;
log_angle( log_angle > d ) = log_angle( log_angle > d ) - d;

rn = [ log2( 0.5 ), log2( R1 ) ];
log_amp( log_amp <= rn( 1 ) | log_amp > rn( 2 ) ) = 100;
log_amp( log_amp ~= 100 ) = 1;
log_amp( log_amp == 100 ) = 0;

log_angle1 = log_angle .* log_amp;


cos_t = cos(  - a / 180 * pi );
sin_t = sin(  - a / 180 * pi );

d = ba;r = br;
X =  - R2:1:R2;
Y =  - R2:1:R2;
[ XX, YY ] = meshgrid( X, Y );
c_rot = XX * cos_t - YY * sin_t;
r_rot = XX * sin_t + YY * cos_t;

log_angle = atan2( r_rot, c_rot );
log_angle = log_angle / pi * 180;
log_angle( log_angle < 0 ) = log_angle( log_angle < 0 ) + 360;
log_amp = log2( sqrt( c_rot .^ 2 + r_rot .^ 2 ) );

log_angle = round( log_angle * d / 360 );
log_angle( log_angle <= 0 ) = log_angle( log_angle <= 0 ) + d;
log_angle( log_angle > d ) = log_angle( log_angle > d ) - d;

k = zeros( 1, r );
k( 1 ) = 1 / sum( r * d - d + 1 );
for i = 2:r
    k( i ) = k( i - 1 ) + d / sum( r * d - d + 1 );
end
k = sqrt( k );
rn = [ log2( 0.5 ), log2( R2 * k ) ];
log_amp( log_amp <= rn( 1 ) | log_amp > rn( r + 1 ) ) = 100;
for i = 1:r
    log_amp( log_amp > rn( i ) & log_amp <= rn( i + 1 ) ) = i;
end
log_amp( log_amp == 100 ) = 0;

log_angle2 = log_angle;
log_amp2 = log_amp;
end
