% clear all
clearvars
clear
clc
clearvars -GLOBAL
close all
format shorte

global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      %metres (32.1740 ft) per s²
C.am = 1.66053892e-27;

TotalElectrons = 10000;
T= 300;
Tau = 0.2e-12;
length = 200e-9;
height = 100e-9;
Px= zeros(1,TotalElectrons);
NewPx= zeros(1,TotalElectrons);
Py= zeros(1,TotalElectrons);
NewPy= zeros(1,TotalElectrons);
Vx= zeros(1,TotalElectrons);
Vy= zeros(1,TotalElectrons);

Vth = sqrt((C.kb*T)/(C.m_0*0.26));
time = 1000;
dt = 15e-15;
lambda = Vth*Tau;
fprintf('The thermal velocity is equal to %d\n',Vth)
fprintf('The mean free path is equal to %d\n',lambda)

figure('Name','Electron Paths')

Px = rand(1, TotalElectrons)*200e-9;       %Generate the random x location
Py = rand(1, TotalElectrons)*100e-9;       %generate the random y location

RandAng = rand(1, TotalElectrons)*2*pi;
Vx = Vth*sin(RandAng);
Vy = Vth*cos(RandAng);

hold on

for loops=1:time
    NewPx = Vx*dt+Px;
    NewPy = Vy*dt+Py;
    
    
    ix = NewPx>length;
    NewPx(ix) = NewPx(ix)-length;
    Px(ix) = Px(ix)-length;
    
    ix = NewPx<0;
    NewPx(ix) = NewPx(ix)+length;
    Px(ix) = Px(ix)+length;

    ix = NewPy<0;
    Vy(ix) = -Vy(ix);
    
    ix = NewPy>height;
    Vy(ix) = -Vy(ix);
    
    subplot (2,1,1)
    plot([Px(1) NewPx(1)], [Py(1) NewPy(1)], 'b')
    plot([Px(2) NewPx(2)], [Py(2) NewPy(2)], 'g')
    plot([Px(3) NewPx(3)], [Py(3) NewPy(3)], 'r')
    plot([Px(4) NewPx(4)], [Py(4) NewPy(4)], 'c')
    plot([Px(5) NewPx(5)], [Py(5) NewPy(5)], 'm')
    plot([Px(6) NewPx(6)], [Py(6) NewPy(6)], 'y')
    hold on
    title('Electron Paths');
    xlim([0 200e-9]);
    ylim ([0 100e-9]);
    pause(0.000001)
    
    Px=NewPx;
    Py=NewPy;

    subplot(2,1,2);
    V = sqrt((Vx.^2)+(Vy.^2));
    Vm = mean(V);
    TCalc = ((Vm^2)*(C.m_0*0.26))/C.kb;
    plot (loops, TCalc, 'b.')
    title ('Temperature (K)');
    xlim([0 time]);
    hold on
end

