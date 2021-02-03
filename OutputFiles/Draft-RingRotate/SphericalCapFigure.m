% from https://www.mathworks.com/matlabcentral/answers/313786-plotting-a-spherical-segment

clear; syms t x y z  real; 
r = 1
a = .043*r    %  40 pct. of r 
b = .130*r    %  60 pct. of r
c = .999*r

eq = x^2+y^2+z^2 == r^2     %  equation of a sphere with center at the origin 
% define high by projecting on xz-plane by setting y=0.  And then x=0  e,g.  the high on  z-axe 
za = solve( subs( eq, [x^2, y^2,r^2], [0,0, (r^2-a^2)] ), z  )     % by Pythagoras, see figure under 
za = za(2)   %  only positive value, upper half of the sphere. 
% Similarly
zb = solve( subs( eq, [x^2, y^2,r^2], [0,0, (r^2-b^2)] ), z  )  % by Pythagoras 
zb =  zb(2)  %  only positive value 
h = za - zb

zc = solve( subs( eq, [x^2, y^2,r^2], [0,0, (r^2-c^2)] ), z  )  % by Pythagoras 
zc =  zc(2)  %  only positive value 

%% plotting 
clear x;
x(t) = sin(t) *r ;      %     -1<x<1 
[X,Y,Z] = sphere;
%s = mesh(X,Z,Y, 'edgecolor', '#A0A0A0', 'FaceAlpha',0.9) ;
s = mesh(X,Y,Z, 'edgecolor', '#C0C0C0', 'FaceAlpha',0.8) ;
s.FaceColor=[1 1 1];
hold on
yApos = sqrt(  r^2 -x^2-za^2); 
yAneg = -yApos;
zA(t) = 0*t + za;
fplot3( x, yApos, zA, [-pi, pi], 'red', 'LineWidth', 3, 'MeshDensity', 720 ), 
fplot3( x, yAneg, zA, [-pi, pi], 'red', 'LineWidth', 3, 'MeshDensity', 720  )
yBpos = sqrt(  r^2 -x^2-zb^2);
yBneg = -yBpos;
zB(t) = 0*t+ zb;
fplot3( x,  yBpos, zB,[-pi, pi], 'red', 'LineWidth', 3, 'MeshDensity', 720  )
fplot3( x,  yBneg, zB, [-pi, pi], 'red', 'LineWidth', 3, 'MeshDensity', 720  )

yCpos = sqrt(  r^2 -x^2-zc^2);
yCneg = -yCpos;
zC(t) = 0*t+ zc;
fplot3( x, yCpos, zC, [-pi, pi*2], 'blue', 'LineWidth', 3, 'MeshDensity', 720  )
fplot3( x, yCneg,  zC, [-pi, pi*2], 'blue', 'LineWidth', 3, 'MeshDensity', 720  )

fplot3( x, yCpos, -zC, [-pi, pi*2], 'blue', 'LineWidth', 3, 'MeshDensity', 720  )
fplot3( x, yCneg, -zC, [-pi, pi*2], 'blue', 'LineWidth', 3, 'MeshDensity', 720  )


xlabel('x'), ylabel('y'), zlabel('z'), axis equal , grid off 
axis( [-1,1, -1,1,  -1,1] *r ) 
campos( [5,5,5])