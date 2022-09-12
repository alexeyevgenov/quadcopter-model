

global Jr Ix Iy Iz b d l m g Kpz Kdz Kpp Kdp Kpt Kdt Kpps Kdps ZdF PhidF ThetadF PsidF ztime phitime thetatime psitime Zinit Phiinit Thetainit Psiinit Uone Utwo Uthree Ufour Ez Ep Et Eps

%koefficienti inercii regulyatorov dlya 4 regima%
%tangaj, kren, riskanie, navesanie%


kpp = 0.8;   
kdp = 0.4;

kpt = 1.2;
kdt = 0.4;

kpps = 1;
kdps = 0.4;

kpz = 100;
kdz = 80;

Gains = [kpp kdp kpt kdt kpps kdps kpz kdz];

disp(Gains);

% Quadrotor constants
Ix = 7.5*10^(-3);  % moment inercii kvadrokopterd/ proeksiya x
Iy = 7.5*10^(-3);  % moment inercii kvadrokopterd/ proeksiya y
Iz = 1.3*10^(-2);  % moment inercii kvadrokopterd/ proeksiya z
Jr = 6.5*10^(-5)+0.1;  % moment inercii vokrug oceipropellera 
b = 3.13*10^(-5);  % Thrust factor
d = 7.5*10^(-7);  % Drag factor
l = 0.5;  % rasstoyanie do centra kvadrokoptera
m = 0.65;  % Massa
g = 9.81;   % uskorenie svobodnogo padeniya

% % Controlling the Quadrotor
sim('quadsim')
