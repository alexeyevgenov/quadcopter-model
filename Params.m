% Quadrotor constants

%global quad.b quad.k quad.l quad.m quad.g
%-----------------------------------
quad.ro=1.184; % air density

quad.m=1.284;%1.282;
quad.g=9.81;

quad.b=0.17246e-5;
quad.k=0.6633e-7;

%quad.b=0.9937e-8;
%quad.k=quad.b/26;

quad.l=0.254;

quad.Ixx = 8.283e-2;
quad.Iyy = 8.283e-2;
quad.Izz = 8.283e-2;
quad.J = [quad.Ixx 0 0; 
    0 quad.Iyy 0; 
    0 0 quad.Izz];


%quad.Cp=quad.Ct*sqrt(quad.Ct/2);
%quad.S=quad.r^2*pi;
%quad.b=quad.Ct*quad.ro*quad.S*quad.r.^2;
%quad.k=quad.Cp*quad.ro*quad.S*quad.r.^3;

%-----------------------

quad.verbose = false;
%sim('IRIS')