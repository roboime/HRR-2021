if 0 %0 para rodar alg aprendizado
    clear all
    close all
    clc
    robosimulador_DataFile;
    Input_Simulink;
    gait=0.6;
    gait_c=1.4;
    passos=4;
    ang_tr=16;
    ang_q=19;
    axf=80;
    azf=80;
end
%taxa=1/400; %Taxa de move 
delta=0.0001;
graf=0;
grafs=0;
csv=0;
for cont_passos=1:passos
    for j=1:round(gait_c/taxa) 
        Qc(round(gait_c/taxa)*(cont_passos-1)+j+1,1)=-ang_q*sin(pi*taxa*j/gait_c)*(-1)^(1-mod(cont_passos,2));
        Qc(round(gait_c/taxa)*(cont_passos-1)+j+1,2)=ang_q*sin(pi*taxa*j/gait_c)*(-1)^(1-mod(cont_passos,2));
        Qc(round(gait_c/taxa)*(cont_passos-1)+j+1,3)=-ang_q*sin(pi*taxa*j/gait_c)*(-1)^(1-mod(cont_passos,2));
        Qc(round(gait_c/taxa)*(cont_passos-1)+j+1,4)=ang_q*sin(pi*taxa*j/gait_c)*(-1)^(1-mod(cont_passos,2));
    end
end
for cont_passos=1:passos
    for j=1:round(gait/taxa) 
        if cont_passos==1
            Q_tr(round(gait/taxa)*(cont_passos-1)+j+1,1)=-ang_tr*sin(0.5*pi*taxa*j/gait)*(-1)^(1-mod(cont_passos,2));
            Q_tr(round(gait/taxa)*(cont_passos-1)+j+1,2)=ang_tr*sin(0.5*pi*taxa*j/gait)*(-1)^(1-mod(cont_passos,2));
        else
            Q_tr(round(gait/taxa)*(cont_passos-1)+j+1,1)=-ang_tr*sin(-pi/2+pi*taxa*j/gait)*(-1)^(1-mod(cont_passos,2));
            Q_tr(round(gait/taxa)*(cont_passos-1)+j+1,2)=ang_tr*sin(-pi/2+pi*taxa*j/gait)*(-1)^(1-mod(cont_passos,2));
        end
    end
end
Cinematica_Heli_Sagital;
for cont_passos=1:passos
    for j=1:round(gait_c/taxa)/2 
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),1)=Qc(j+round(gait_c/taxa)*(cont_passos-1),1);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),5)=Qc(j+round(gait_c/taxa)*(cont_passos-1),2);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),8)=Qc(j+round(gait_c/taxa)*(cont_passos-1),3);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),12)=Qc(j+round(gait_c/taxa)*(cont_passos-1),4);
        if cont_passos~=1
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),2)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),2);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),3)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),3);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),4)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),4);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),6)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),6);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),7)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),7);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),9)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),9);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),10)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),10);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),11)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),11);
        end
    end
    for j=(round(gait_c/taxa)/2+1):(round(gait_c/taxa)/2+round(gait/taxa))
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),2)=Q(j-round(gait_c/taxa)/2+round(gait/taxa)*(cont_passos-1),1);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),3)=Q(j-round(gait_c/taxa)/2+round(gait/taxa)*(cont_passos-1),2);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),4)=Q(j-round(gait_c/taxa)/2+round(gait/taxa)*(cont_passos-1),3);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),9)=Q(j-round(gait_c/taxa)/2+round(gait/taxa)*(cont_passos-1),4);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),10)=Q(j-round(gait_c/taxa)/2+round(gait/taxa)*(cont_passos-1),5);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),11)=Q(j-round(gait_c/taxa)/2+round(gait/taxa)*(cont_passos-1),6);
        
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),6)=Q_tr(j-round(gait_c/taxa)/2+round(gait/taxa)*(cont_passos-1),1);
        Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),7)=Q_tr(j-round(gait_c/taxa)/2+round(gait/taxa)*(cont_passos-1),2);
        
          if 1
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),1)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),1);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),5)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),5);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),8)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),8);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),12)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),12);
        end
    end
    for j=(round(gait_c/taxa)/2+round(gait/taxa)+1):(round(gait_c/taxa)+round(gait/taxa)) 
            if mod(cont_passos,2)==1
                Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),1)=Qc(j-round(gait/taxa),1);
                Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),5)=Qc(j-round(gait/taxa),2);
                Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),8)=Qc(j-round(gait/taxa),3);
                Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),12)=Qc(j-round(gait/taxa),4);
            end
            if mod(cont_passos,2)==0
                Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),1)=-Qc(j-round(gait/taxa),1);
                Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),5)=-Qc(j-round(gait/taxa),2);
                Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),8)=-Qc(j-round(gait/taxa),3);
                Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),12)=-Qc(j-round(gait/taxa),4);
            end
        %Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),1)=Qc(j-round(gait/taxa),1);
        %Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),5)=Qc(j-round(gait/taxa),2);
        %Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),8)=Qc(j-round(gait/taxa),3);
        %Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),12)=Qc(j-round(gait/taxa),4);
        if 1
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),2)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),2);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),3)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),3);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),4)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),4);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),6)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),6);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),7)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),7);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),9)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),9);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),10)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),10);
            Qr(j+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),11)=Qr(j-1+(round(gait/taxa)+round(gait_c/taxa))*(cont_passos-1),11);
        end
    end
    
end
clear  Q X_E Z_E X_D Z_D X Y Z Tr 
Q=Qr;
%% Cinem�tica (Desenvoler o movimento) 
lx=[-64    0  21    0  0 0 0     0 0 0   -21    0 64    ]; %Comprimentos na direção x, ultimo � transforma��o da ferramenta
ly=[0      0  0     0  0 0 82.09 0 0 0   0      0 0     ]; %Comprimentos na dire��o y
lz=[19.37  0  70.38 82 0 0 0     0 0 -82 -70.38 0 -19.37]; %Comprimentos na dire��o z
lxr(1)=lx(1);
lyr(1)=ly(1);
lzr(1)=lz(1);
for i=2:length(lz) %Resultante dos comprimentos para modelagem por helicoide
    lxr(i)=lx(i)+lxr(i-1);
    lyr(i)=ly(i)+lyr(i-1); 
    lzr(i)=lz(i)+lzr(i-1);
end
q=[0 0 0 0 0 0 0 0 0 0 0 0]; %angulos
X_E(1)=0;
Z_E(1)=0;
X_D(1)=0;
Z_D(1)=0;
    X=zeros(1,14);
    Y=zeros(1,14);
    Z=zeros(1,14);
    X(1)=0;
    Y(1)=0;
    Z(1)=0;
    X(2)=lxr(1);
    Y(2)=lyr(1);
    Z(2)=lzr(1);
    fx=1;%rand(1);
    fz=1;%rand(1);
    for cont_passos=1:passos
        if 1
            q_aux=q;
            for i=1:length(q_aux)
                q(i)=-q_aux(length(q_aux)+1-i);                
            end
        end
        q_aux=q;
        if cont_passos>0
            lx=[X(14)-64     0  21    0  0 0 0                          0 0 0   -21    0 64    ]; %Comprimentos na direção x, ultimo � transforma��o da ferramenta
            ly=[Y(14)        0  0     0  0 0 82.09*(-1)^(cont_passos-1) 0 0 0   0      0 0     ]; %Comprimentos na dire��o y
            lz=[Z(14)+19.37  0  70.38 82 0 0 0                          0 0 -82 -70.38 0 -19.37]; %Comprimentos na dire��o z
            lxr(1)=lx(1);
            lyr(1)=ly(1);
            lzr(1)=lz(1);
            X(2)=lxr(1);
            Y(2)=lyr(1);
            Z(2)=lzr(1);
            for i=2:length(lz) %Resultante dos comprimentos para modelagem por helicoide
                lxr(i)=lx(i)+lxr(i-1);
                lyr(i)=ly(i)+lyr(i-1);
                lzr(i)=lz(i)+lzr(i-1);
            end
            PH=[-1^(cont_passos-1) 0                  0                   lxr(1)  lyr(1)  lzr(1) %sx sy sz s0x s0y s0z(Parametros helicoide)
                0                -1^(cont_passos-1)	0                   lxr(2)  lyr(2)  lzr(2) %Eixo dos motores de cada junta
                0                -1^(cont_passos-1)	0                   lxr(3)  lyr(3)  lzr(3) 
                0                -1^(cont_passos-1)	0                   lxr(4)  lyr(4)  lzr(4) 
                -1^(cont_passos-1) 0                  0                   lxr(5)  lyr(5)  lzr(5)
                0                0                  -1^(cont_passos-1)    lxr(6)  lyr(6)  lzr(6)
                0                0                  -1^(cont_passos-1)    lxr(7)  lyr(7)  lzr(7)
                -1^(cont_passos-1) 0                  0                   lxr(8)  lyr(8)  lzr(8)
                0                -1^(cont_passos-1)	0                   lxr(9)  lyr(9)  lzr(9) 
                0                -1^(cont_passos-1)	0                   lxr(10)  lyr(10)  lzr(10) 
                0                -1^(cont_passos-1)	0                   lxr(11)  lyr(11)  lzr(11) 
                -1^(cont_passos-1) 0                  0                   lxr(12)  lyr(12)  lzr(12)];
            X(1)=X(14);
            Y(1)=Y(14);
            Z(1)=Z(14);
            X(2)=lxr(1);
            Y(2)=lyr(1);
            Z(2)=lzr(1);
        end
        for j=1:round((gait+gait_c)/taxa) 
            q_aux=q;
            X_aux=X;
            Y_aux=Y;
            Z_aux=Z;
            Tr=eye(4); 
            q=q*pi/180;
            for i=1:length(lz)-1 %Determina a posi��o do rob� (em cada instante de tempo)
                a11=cos(q(i))+(1-cos(q(i)))*PH(i,1)^2;              %Termos da matriz de rota��o da Transformada
                a12=PH(i,1)*PH(i,2)*(1-cos(q(i)))-PH(i,3)*sin(q(i));%por helicoide
                a13=PH(i,1)*PH(i,3)*(1-cos(q(i)))+PH(i,2)*sin(q(i));
                a21=PH(i,2)*PH(i,1)*(1-cos(q(i)))+PH(i,3)*sin(q(i));
                a22=cos(q(i))+(1-cos(q(i)))*PH(i,2)^2;
                a23=PH(i,2)*PH(i,3)*(1-cos(q(i)))-PH(i,1)*sin(q(i));
                a31=PH(i,3)*PH(i,1)*(1-cos(q(i)))-PH(i,2)*sin(q(i));
                a32=PH(i,3)*PH(i,2)*(1-cos(q(i)))+PH(i,1)*sin(q(i));
                a33=cos(q(i))+(1-cos(q(i)))*PH(i,3)^2;
                T=[a11 a12 a13 PH(i,4)-PH(i,4)*a11-PH(i,5)*a12-PH(i,6)*a13 %Matriz Tab
                    a21 a22 a23 PH(i,5)-PH(i,4)*a21-PH(i,5)*a22-PH(i,6)*a23
                    a31 a32 a33 PH(i,6)-PH(i,4)*a31-PH(i,5)*a32-PH(i,6)*a33
                    0 0 0 1];
                Tf=[1 0 0 lxr(i+1)   %Matriz de transforma��o da ferramenta(prox junta)
                    0 1 0 lyr(i+1)
                    0 0 1 lzr(i+1) 
                    0 0 0 1];
                Tr=Tr*T;
                TF=Tr*Tf; %Transforma�ao final(da prox junta)
                X(i+2)=TF(1,4);  %posi��o x
                Y(i+2)=TF(2,4); %posi��o y
                Z(i+2)=TF(3,4); %posi��o z
                Theta_X(i+1)=atan(TF(3,2)/TF(3,3));
                Theta_Y(i+1)=atan(-TF(3,1)/sqrt(TF(3,2)^2+TF(3,3)^2));
                Theta_Z(i+1)=atan(TF(2,1)/TF(1,1));
            end 
            Xg(j)=X(i+2);
            Yg(j)=Y(i+2);
            Zg(j)=Z(i+2);
            if  grafs && (mod(j,20)==0 ||j==1)
               
                plot3(X,Y,Z);
                hold on
                grid on
                grid minor
                axis([-100 300 -100 300 -100 300])
                %M(2*i)=getframe;
                pause(taxa*10)
                %cla;
                hold off
            end    
            for k=1:12
                q(k)=Q(j+(round((gait+gait_c)/taxa))*(cont_passos-1),k);
            end
          
      
            if mod(cont_passos,2)==1
                X_E(round((gait+gait_c)/taxa)*(cont_passos-1)+j+1)=X(length(X));
                Z_E(round((gait+gait_c)/taxa)*(cont_passos-1)+j+1)=Z(length(Z));
                X_D(round((gait+gait_c)/taxa)*(cont_passos-1)+j+1)=X_D(round((gait+gait_c)/taxa)*(cont_passos-1)+j);
                Z_D(round((gait+gait_c)/taxa)*(cont_passos-1)+j+1)=Z_D(round((gait+gait_c)/taxa)*(cont_passos-1)+j);
            end
            if mod(cont_passos,2)==0
                X_E(round((gait+gait_c)/taxa)*(cont_passos-1)+j+1)=X_E(round((gait+gait_c)/taxa)*(cont_passos-1)+j);
                Z_E(round((gait+gait_c)/taxa)*(cont_passos-1)+j+1)=X_E(round((gait+gait_c)/taxa)*(cont_passos-1)+j);
                X_D(round((gait+gait_c)/taxa)*(cont_passos-1)+j+1)=X(length(X));
                Z_D(round((gait+gait_c)/taxa)*(cont_passos-1)+j+1)=Z(length(Z));
            end
        end
    end
%Q_aux=Q;
%for i=1:length(Q_aux)
    %Q(i+length(Q_aux),8:-1:1)=Q_aux(i,1:8);
%end
if 1
Q_aux=Q;
tam=length(Q);
for i=1:round(gait/taxa)
    for j=1:12
        Q(tam+i,j)=Q_aux(tam,j)-i*Q_aux(tam,j)/round(gait/taxa);
    end
end
end
if  1 %Ativa espelhamento de passo
    Q_aux=Q;
    tam=length(Q);
    for i=1:tam
        for j=1:12
            Q(tam+i,j)=Q_aux(i,13-j);
        end
    end
end
if 0 %Ativa espelhamento de passo - por cima
    Q_aux=Q;
    tam=length(Q);
    for i=1:tam
        for j=1:12
            Q(i,j)=Q_aux(i,13-j);
        end
    end
end
if graf
    figure
    plot(X_E,Z_E)
    hold on
    plot(X_D,Z_D) 
    figure
    for j=1:length(Q)
        Q_TD_E(j)=Q(j,1);
        Q_TD_F(j)=Q(j,2);
        Q_JD(j)=Q(j,3);
        Q_QD_F(j)=Q(j,4);
        Q_QD_E(j)=Q(j,5);
        Q_QD_R(j)=Q(j,6);
        Q_QE_R(j)=Q(j,7);
        Q_QE_E(j)=Q(j,8);
        Q_QE_F(j)=Q(j,9);
        Q_JE(j)=Q(j,10);
        Q_TE_F(j)=Q(j,11);
        Q_TE_E(j)=Q(j,12);
        Q_TIME(j)=j*taxa;
    end
    plot(Q_TIME,Q_TD_E);
    hold on
    plot(Q_TIME,Q_TD_F);
    plot(Q_TIME,Q_JD);
    plot(Q_TIME,Q_QD_F);
    %plot(Q_TIME,Q_QD_E);
    plot(Q_TIME,Q_QD_R);
    plot(Q_TIME,Q_QE_R);
    %plot(Q_TIME,Q_QE_E);
    plot(Q_TIME,Q_QE_F);
    plot(Q_TIME,Q_JE);
    plot(Q_TIME,Q_TE_F);
    plot(Q_TIME,Q_TE_E);
    %grid on
    %grid minor    
    legend('Tonozelo D Eversão - Quadril E Eversão','Tonozelo D Flexão','Joelho D','Quadril D flexão','Quadril D Rotação','Quadril E Rotação','Quadril E Flexão','Joelho E','Tornozelo E Flexão','Tornozelo E eversão - Quadril eversão D');
    xlabel('Seconds');
    ylabel('Degrees');
    RPM(1,1:12)=0;
    for i=2:length(Q)
        for j=1:12
            RPM(i,j)=(Q(i,j)-Q(i-1,j))*60/(taxa*360);
            if RPM(i,j)>50
                RPM(i,j)=50;
            end
            if RPM(i,j)<-50
                RPM(i,j)=-50;
            end
            if RPM(i,j)<0
                RPM(i,j)=-RPM(i,j);
            end
        end
    end
    figure
    plot(RPM)
    M=[0 0];
    if 1 %Completo
         for i=1:length(Q)
            for j=1:12
                M((i-1)*12+j,1)=Q(i,j);
                M((i-1)*12+j,2)=RPM(i,j);
                %fprintf("%f,%f\n",Q(i,j),RPM(i,j));
            end
        end
    end
    if 0 %inicio
        for i=1:4*(length(Q)-1)/17
            for j=1:4
                M((i-1)*4+j,1)=Q(i,j);
                M((i-1)*4+j,2)=RPM(i,j);
                %fprintf("%f,%f\n",Q(i,j),RPM(i,j));
            end
        end
    end
    if 0%meio
        for i=1:4*length(Q)/17
            for j=1:6
                M((i-1)*6+j,1)=Q(i+4*(length(Q)-1)/17,j);
                M((i-1)*6+j,2)=RPM(i+4*(length(Q)-1)/17,j);
                %fprintf("%f,%f\n",Q(i,j),RPM(i,j));
            end
        end
    end
    if 0 %fim
        for i=1:4*length(Q)/17
            for j=1:6
                %M((i-1)*6+j,1)=Q(i+13*(length(Q)-1)/17,j);
                %M((i-1)*6+j,2)=RPM(i+13*(length(Q)-1)/17,j);
                %fprintf("%f,%f\n",Q(i,j),RPM(i,j));
            end
        end
    end
    if csv
        strcsv=strcat('V4_200Hz_','G',num2str(gait),'-X',num2str(axf),'-Z',num2str(azf),'-Tr',num2str(ang_tr)','-P12','-Pc0.1','-T-1','.csv');
        %csvwrite(strcsv,M);
        csvwrite('Passo_1_G0.6_Gc1.4_Ax80_Az80_AngTr16_AngQ19.csv',M);
    end    
end
j=length(Q); %ÚLTIMAS LINHAS DO CÓDIGO
for i=0:j-1 
    Time(i+1)=i*taxa;
end 