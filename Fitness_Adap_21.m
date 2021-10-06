%% Fitness Adaption
minF=10e5;
maxF=-10e5;
for cont_p=1:length(Pp) 
    clear Q Qr Qc Q_tr Time
    gait=Pp(cont_p,1)*gait_B;
    gait_c=Pp(cont_p,1)*gait_c_B;
    ang_tr=Pp(cont_p,2);
    ang_q=Pp(cont_p,3);
    axf=Pp(cont_p,4); 
    azf=Pp(cont_p,5);  
    
    
    if (cont==1 && cont_p==1) % se eh o primeiro elemento da primeira geracao
        taxa=1/400; % frequencia de envio de comandos motores
        Cinematica_Comb; 
        simout = sim('robosimulador','FastRestart','off','ReturnWorkspaceOutputs','on');
        TP=(gait+gait_c)/taxa; %fator de tamanho de matriz
    else 
        taxa=(gait+gait_c)/TP;
        Cinematica_Comb; 
        simout = sim('robosimulador','FastRestart','on','ReturnWorkspaceOutputs','on');
    end 
    
    %Funcao de maximacao
    dist_x=-1; % garantir o minimo
    queda=simout.Queda(length(simout.Queda));
    for cont_simout=1:length(simout.X)
        if dist_x<simout.X(cont_simout)
            dist_x=simout.X(cont_simout);
            p_xmax=cont_simout;
        end
    end
    %pegando as variavei exportadas pelo simulink
    dist_y=simout.Y(p_xmax);
    tempo=simout.Tempo(length(simout.Tempo));
    t1max=simout.T1max(length(simout.T1max));
    t2max=simout.T2max(length(simout.T2max));
    t3max=simout.T3max(length(simout.T3max));
    t4max=simout.T4max(length(simout.T4max));
    t5max=simout.T5max(length(simout.T5max));
    t6max=simout.T6max(length(simout.T6max));
    t7max=simout.T7max(length(simout.T7max));
    t8max=simout.T8max(length(simout.T8max));
    %Primeira parte da penalização desvio segundo leva em consideração
    %queda e torque
    F=0;
    %  distância  velocidade       (penalidade desvio )               (penalidade queda e torque alto)    
    F=20*dist_x+12*dist_x/tempo-2*abs(dist_y/dist_x)-5*(queda+t1max+t2max+t3max+t4max+t5max+t6max+t7max+t8max); %Função recompensa
    
    
    %zerando funcoes de recompenca
    F_V(cont_p)=F;
    if F<minF
       minF=F;
       pmin=cont_p;
    end
    if F>maxF %salvando melhor passo
       maxF=F;
       pmax=cont_p;
       dist_x_m=dist_x;
       dist_y_m=dist_y;
       tempo_m=tempo;
    end
end

Pr=csvread('Parametros.csv');
sizep=size(Pr);
Pr(sizep(1)+1,1)=dist_x_m;
Pr(sizep(1)+1,2)=maxF;
Pr(sizep(1)+1,3)=dist_y_m; %desvio em y (não é angulo)
Pr(sizep(1)+1,4)=Pp(pmax,1)*gait_B;
Pr(sizep(1)+1,5)=Pp(pmax,1)*gait_c_B;
Pr(sizep(1)+1,6)=Pp(pmax,2);
Pr(sizep(1)+1,7)=Pp(pmax,3);
Pr(sizep(1)+1,8)=Pp(pmax,4);
Pr(sizep(1)+1,9)=Pp(pmax,5);
Pr(sizep(1)+1,10)=taxa;
csvwrite('Parametros.csv',Pr); 

%Analisa muta��o
if 1 %ativando a proteço contra a mutação (garante que a mutaço é benéfica (pode covergir em um máximo local))
    if cont>2
        for cont_p=1:length(Pp)
            if F_V_A(cont_p)>F_V(cont_p) 
                Pp(cont_p,1:5)=P_A(cont_p,1:5); 
            end
        end
    end 
end

%Sele��o Natural && Cruzamento
Sum=0;
for cont_p=1:length(Pp)      %Coloca Fmin = 0, adicionando no restante
    if (cont_p ~=pmax)
        F_V_s(cont_p)=F_V(cont_p)-minF;
        Sum=Sum+F_V_s(cont_p);
    end
end
for cont_p=1:length(Pp)      %Calcula propabilidade de escolha para elementos restantes
    if (cont_p ~=pmax)
        P_F(cont_p)=F_V_s(cont_p)/Sum;
    end    
end
R=round(1);
roll=0;
for cont_p=1:10   %Escolhe parceiro
    if (cont_p ~=pmax)
        R=R-P_F(cont_p);
        roll=roll+1;
    end    
    if R<0
       break; %[] trocar por while
    end
end


Pp(pmin,1:5)=(Pp(roll,1:5)+Pp(pmax,1:5))/2; %cruzamento melhor com o ganhador da roleta

%Pp(pmin,1:3)=Pp(roll,1:3);
 

