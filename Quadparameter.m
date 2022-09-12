% Quadrotor constants

%global quad.b quad.k quad.l quad.m quad.g
%-----------------------------------
quad.ro=1.184; % air density
quad.r=0.165;%0.127;  %prop diameter
quad.Ct=0.0048; %0.001396; 

quad.m=0.468;%1.282;
quad.g=9.81;

quad.b=2.98e-6;
quad.k=1.14e-7;

quad.l=0.225;

quad.Ixx = 4.856e-3;
quad.Iyy = 4.856e-3;
quad.Izz = 8.801e-3;
quad.J = [quad.Ixx 0 0; 
    0 quad.Iyy 0; 
    0 0 quad.Izz];

Ax = 0.25;
Ay = 0.25;
Az = 0.25;
quad.A = [Ax 0 0; 
    0 Ay 0; 
    0 0 Az];

quad.Jr=3.375e-5;

%quad.Cp=quad.Ct*sqrt(quad.Ct/2);
%quad.S=quad.r^2*pi;
%quad.b=quad.Ct*quad.ro*quad.S*quad.r.^2;
%quad.k=quad.Cp*quad.ro*quad.S*quad.r.^3;

%-----------------------

quad.verbose = false;
%sim('Allgemeine_Schematische')