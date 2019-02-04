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
C.meff = C.m_0 *0.26;
C.am = 1.66053892e-27;

TotalElectrons = 10000;


Tau = 0.2e-12;  
length = 200e-9;
height = 100e-9;
time = 500;
dt = 15e-15;
OGTemp = 300;

%Create Vectors Sized Based on Electrons
Px= zeros(1,TotalElectrons);
NewPx= zeros(1,TotalElectrons);
Py= zeros(1,TotalElectrons);
NewPy= zeros(1,TotalElectrons);
Vx= zeros(1,TotalElectrons);
Vy= zeros(1,TotalElectrons);
Scat = zeros(1,TotalElectrons);
PScat = (1-exp(-dt/Tau));
TotalTemp = 0;
AverageTemp = 0;

%Calculate the initial thermal voltage
Vth = sqrt(2*(C.kb*300)/(C.m_0*0.26));
lambda = Vth*Tau;

%Output the values for initial thermal velocity and MFP.
fprintf('The initial thermal velocity is equal to %d\n',Vth)
fprintf('The mean free path is equal to %d\n',lambda)

%Name the figure window
figure('Name','Electron Paths')

%Calculate the random starting position for every electron
Px = rand(1, TotalElectrons)*200e-9;       %Generate the random x location
Py = rand(1, TotalElectrons)*100e-9;       %generate the random y location

%Pick a random velocity in X and Y from the distribution
Vx = (sqrt(C.kb*300/(C.m_0*0.26))*randn(1,TotalElectrons)); %Generate random x velocity
Vy = (sqrt(C.kb*300/(C.m_0*0.26))*randn(1,TotalElectrons)); %Generate random y velocity

%Calculate the Voltage
V2 = (Vx.^2)+(Vy.^2);
V = sqrt(V2);


subplot(2,2,3)
histogram(V,40);                        %Plot histogram of velocity
title ('Distribution of Random Velocities');
Vm = mean(V2);                          %Find the average velocity of electrons
TCalc = ((Vm)*(C.m_0*0.26))/2/C.kb;     %Calculate the temperature


for loops=1:time
    Scat = rand(1,TotalElectrons);
    
    Vx(Scat<PScat) = (Vth/sqrt(2))*randn(1);
    Vy(Scat<PScat) = (Vth/sqrt(2))*randn(1);

    NewPx = Vx*dt+Px;                       %Calculate new x position
    NewPy = Vy*dt+Py;                       %Calculate new y position
    
    %Check right boundary
    ix = NewPx>length;
    NewPx(ix) = NewPx(ix)-length;
    Px(ix) = Px(ix)-length;
        
    %Check left boundary
    ix = NewPx<0;
    NewPx(ix) = NewPx(ix)+length;
    Px(ix) = Px(ix)+length;

    %Check bottom boundary
    ix = NewPy<0;
    Vy(ix) = -Vy(ix);
    
    %Check top boundary
    ix = NewPy>height;
    Vy(ix) = -Vy(ix);
    
    %Plot electrons
    subplot (2,2,1)
    
    plot([Px(1) NewPx(1)], [Py(1) NewPy(1)], 'b')       %Electron 1
    plot([Px(2) NewPx(2)], [Py(2) NewPy(2)], 'g')       %Electron 2
    plot([Px(3) NewPx(3)], [Py(3) NewPy(3)], 'r')       %Electron 3 
    plot([Px(4) NewPx(4)], [Py(4) NewPy(4)], 'c')       %Electron 4
    plot([Px(5) NewPx(5)], [Py(5) NewPy(5)], 'm')       %Electron 5
    plot([Px(6) NewPx(6)], [Py(6) NewPy(6)], 'y')       %Electron 6
    plot([Px(7) NewPx(7)], [Py(7) NewPy(7)], 'k')       %Electron 6
    
    hold on
    title('Electron Paths');
    xlim([0 200e-9]);
    ylim ([0 100e-9]);
    pause(0.000001)
    
    %Update old positions
    Px=NewPx;
    Py=NewPy;

    %Calculate and plot temperature
    V2 = (Vx.^2)+(Vy.^2);
    V = sqrt(V2);
         
    Vm = mean(V2);
    TCalc = ((Vm)*(C.m_0*0.26))/2/C.kb;
    subplot(2,2,2)
    plot([loops-1 loops], [OGTemp TCalc],'b');
    TotalTemp = TotalTemp+TCalc;
    AverageTemp = (TotalTemp/loops);
    
    legend (strcat('Temp; Average Temp: ', num2str(round(AverageTemp)),' k'));
    title ('Temperature (K)');
    xlim([0 time]);
    ylim ([0 600]);
    OGTemp = TCalc;
    hold on
end

