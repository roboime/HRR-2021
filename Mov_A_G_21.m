%% Inicializacao
clc; clear all;
robosimulador_DataFile; %Robot's Geometry
Input_Simulink; %Enviroment parameters 

%% Inputs 
passos = 2; %Numero de passos
gait_B=0.6; % movimento longitudinal (segundos)
gait_c_B=1.4; % movimento de roll (segundos)

%Pp - População inicial         [] implementar busca de arquivo Pp
% K de tempo | Angulo de rot torso | Angulo de Inclinacao | x elicse | z elipse
Pp=[1.0 16 19 80 80
    0.7 10 12 60 30
    0.5 12 14 60 50
    2.0 10 15 100 60
    0.9 11 20 99 60
    0.7 8  14 75 40
    0.9 9  12 62 55
    1.2 6  13 65 65
    1.0 20 11 50 70
    0.9 22 19 44 74];

%% Processo
for  cont=1:150
    % [] Adicionar condicao de parada
    Fitness_Adap_21;
    mutation_21;
    cont
end 
 