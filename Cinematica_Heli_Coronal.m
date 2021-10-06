%HELICOIDES
clear all
clc
close all
ayf=0;
azf=100;
%axt=20;
%azt=20;
Pc_g=0;
gait=0.4;
%ang_tr=8;
passos=1;
graf=1; % quanto quer mostrar a cinemática graf = 1
csv=0; %csv = 1 para gerar o arquivo de memória pra o robô real
lx=[-64   21       0         -21  	 64   ]; %Comprimentos na direção x, ultimo � transforma��o da ferramenta
ly=[0     0        82.09	 0       0    ]; %Comprimentos na dire��o y
lz=[19.37 152.38   0        -152.38 -19.37]; %Comprimentos na dire��o z
lxr(1)=lx(1);
lyr(1)=ly(1);
lzr(1)=lz(1);
for i=2:length(lz) %Resultante dos comprimentos para modelagem por helicoide
    lxr(i)=lx(i)+lxr(i-1);
    lyr(i)=ly(i)+lyr(i-1); 
    lzr(i)=lz(i)+lzr(i-1);
end
PH=[-1	0	0  lxr(1)  lyr(1)  lzr(1) %sx sy sz s0x s0y s0z(Parametros helicoide)
    1	0	0  lxr(2)  lyr(2)  lzr(2) %Eixo dos motores de cada junta
    1	0	0  lxr(3)  lyr(3)  lzr(3) 
    -1	0	0  lxr(4)  lyr(4)  lzr(4)];
%% Cinem�tica (Desenvoler o movimento) 
q=[0 0 0 0 ]; %angulos 
Q(1,1:4)=q;
taxa=1/400; %Taxa de move 
delta=0.0001;  
delay_ini=0; 
X_E(1)=0;
Z_E(1)=0;
X_D(1)=0;
Z_D(1)=0;
    X=zeros(1,6);
    Y=zeros(1,6);
    Z=zeros(1,6);
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
            lx=[X(6)-64    21       0        -21  	 64   ]; %Comprimentos na direção x, ultimo � transforma��o da ferramenta
            ly=[Y(6)       0        82.09*(-1)^(cont_passos-1)	 0       0    ]; %Comprimentos na dire��o y
            lz=[Z(6)+19.37 152.38   0        -152.38 -19.37]; %Comprimentos na dire��o z
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
            PH=[-1^(cont_passos) 0	0  lxr(1)  lyr(1)  lzr(1) %sx sy sz s0x s0y s0z(Parametros helicoide)
                -1^(cont_passos) 0	0  lxr(2)  lyr(2)  lzr(2) %Eixo dos motores de cada junta
                -1^(cont_passos) 0	0  lxr(3)  lyr(3)  lzr(3) 
                -1^(cont_passos) 0	0  lxr(4)  lyr(4)  lzr(4)];
            X(1)=X(6);
            Y(1)=Y(6);
            Z(1)=Z(6);
            X(2)=lxr(1);
            Y(2)=lyr(1);
            Z(2)=lzr(1);
        end
        for j=1:round(gait/taxa) 
            q_aux=q;
            if j<Pc_g*round(gait/taxa)
                Xd=[0 0 0 0];
                %Xd=[axf*sin(pi*j*taxa/gait) 0 azf*cos(pi*j*taxa/gait)  0 0 0 0 0];  
            else 
                 %Xd=[ayf*cos(pi*(j-Pc_g*(gait/taxa))*taxa*(1-Pc_g)/gait)  azf*cos(pi*(j-Pc_g*(gait/taxa))*taxa*(1-Pc_g)/gait) 0  0 ]; 
                 Xd=[ayf azf*cos(pi*(j-Pc_g*(gait/taxa))*taxa*(1-Pc_g)/gait) 0  0 ]; 
            end            
            
           %Xd=[0.05 0.05 0  0 0 0 ]; %a*2*fCMz*cos(pi*j*taxa/gait Velocidade
            Xd=transpose(Xd);
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
            if  graf && (mod(j,10)==0 ||j==1)
                plot3(X,Y,Z);
                hold on
                grid on
                grid minor
                axis([-100 400 -100 300 -100 350])
                %M(2*i)=getframe;
                pause(taxa*500)
                %cla;
                %hold off
            end 
            for k=1:length(q)
                q(k)=q(k)+delta;
                Tr=eye(4); 
                for i=1:length(lz)-1 %Determina a posi��o do rob�(em cada instante de tempo
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
                    Tf=[1 0 0 lxr(i+1)   %Matriz de transforma��o da ferramenta
                        0 1 0 lyr(i+1)
                        0 0 1 lzr(i+1) 
                        0 0 0 1];
                    Tr=Tr*T;
                    TF=Tr*Tf; %Transforma�ao final
                    Xj(i+1)=TF(1,4);  %posi��o x
                    Yj(i+1)=TF(2,4); %posi��o y
                    Zj(i+1)=TF(3,4); %posi��o z
                    Theta_Xj(i+1)=atan(TF(3,2)/TF(3,3));
                    Theta_Yj(i+1)=atan(-TF(3,1)/sqrt(TF(3,2)^2+TF(3,3)^2));
                    Theta_Zj(i+1)=atan(TF(2,1)/TF(1,1));
                end
                J(1,k)=(Yj(i+1)-Y(i+1))/delta;
                J(2,k)=(Zj(i+1)-Z(i+1))/delta; 
                J(3,k)=(Theta_Xj(i+1)-Theta_X(i+1))/delta;
                J(4,k)=(Theta_Xj(4)-Theta_Y(4))/delta;
                %J(7,k)=(CMX_Fj-CMX_F)/delta; 
                %J(7,k)=(CMY_Fj-CMY_F)/delta;
                %J(5,k)=(Theta_Yj(3)-Theta_Y(3))/delta;
                %J(9,k)=(CMZ_Fj-CMZ_F)/delta;
                %J(6,k)=(Theta_Yj(4)-Theta_Y(4))/delta;
                %J(10,k)=(Theta_Zj(i-5)-Theta_Z(i-5))/delta;                 
                q(k)=q(k)-delta;
            end
            Jq=inv(J); 
            %det(Jq)
            %det(J)
            qd=Jq*Xd;
            qd=transpose(qd);
            q=qd*taxa+q;
            q=q*180/pi;
   
            if j>Pc_g*round(gait/taxa)
                if j>2
                    q(1)=q(1)-q_aux(1)/round(gait/taxa);
                    q(2)=q(2)-q_aux(2)/round(gait/taxa);
                end
            end
            if q(2)<0
                q(2)=-q(2);
            end            
            q(1)=-q(2);
            
            %q(4)=-q(3);
            for cont_ver=1:4              
                if (q(cont_ver)-q_aux(cont_ver))*60/(360*taxa)>50
                    %q(cont_ver)=q_aux(cont_ver)+300*taxa;
                end
                if (q(cont_ver)-q_aux(cont_ver))*60/(360*taxa)<50
                    %q(cont_ver)=q_aux(cont_ver)-300*taxa;
                end
            end
            if mod(cont_passos,2)==1
                for i=1:length(q)
                    if i>2
                        Q(round(gait/taxa)*(cont_passos-1)+j+1,i)=q(i);
                    else
                        Q(round(gait/taxa)*(cont_passos-1)+j+1,i)=-q(i);
                    end
                end
            end
            if mod(cont_passos,2)==0
                for i=1:length(q)
                    if i>2
                        Q(round(gait/taxa)*(cont_passos-1)+j+1,i)=-q(length(q)+1-i);
                    else
                        Q(round(gait/taxa)*(cont_passos-1)+j+1,i)=q(length(q)+1-i);
                    end
                end
            end  
            if mod(cont_passos,2)==1
                X_E(round(gait/taxa)*(cont_passos-1)+j+1)=X(length(X));
                Z_E(round(gait/taxa)*(cont_passos-1)+j+1)=Z(length(Z));
                X_D(round(gait/taxa)*(cont_passos-1)+j+1)=X_D(round(gait/taxa)*(cont_passos-1)+j);
                Z_D(round(gait/taxa)*(cont_passos-1)+j+1)=Z_D(round(gait/taxa)*(cont_passos-1)+j);
            end
            if mod(cont_passos,2)==0
                X_E(round(gait/taxa)*(cont_passos-1)+j+1)=X_E(round(gait/taxa)*(cont_passos-1)+j);
                Z_E(round(gait/taxa)*(cont_passos-1)+j+1)=X_E(round(gait/taxa)*(cont_passos-1)+j);
                X_D(round(gait/taxa)*(cont_passos-1)+j+1)=X(length(X));
                Z_D(round(gait/taxa)*(cont_passos-1)+j+1)=Z(length(Z));
            end
        end
    end
%Q_aux=Q;
%for i=1:length(Q_aux)
    %Q(i+length(Q_aux),8:-1:1)=Q_aux(i,1:8);
%end
Q_aux=Q;
tam=length(Q);
for i=1:round(gait/taxa)
    for j=1:4
        Q(tam+i,j)=Q_aux(tam,j)-i*Q_aux(tam,j)/round(gait/taxa);
    end
end
if 0 %Ativa espelhamento de passo
    Q_aux=Q;
    tam=length(Q);
    for i=1:tam
        for j=1:4
            Q(tam+i,j)=Q_aux(i,9-j);
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
        Q_TD(j)=Q(j,1);
        Q_QD(j)=Q(j,2);
        Q_QE(j)=Q(j,3);
        Q_TE(j)=Q(j,4);
        Q_TIME(j)=j*taxa;
    end
    plot(Q_TIME,Q_TD);
    hold on
    plot(Q_TIME,Q_QD);
    plot(Q_TIME,Q_QE);
    plot(Q_TIME,Q_TE);
    %grid on
    %grid minor    
    legend('Tonozelo D','Quadril eversão D','Quadril eversão E','Tornozelo E');
    xlabel('Seconds');
    ylabel('Degrees');
    RPM(1,1:4)=0;
    for i=2:length(Q)
        for j=1:4
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
    if 0 %Completo
         for i=1:length(Q)
            for j=1:6
                M((i-1)*6+j,1)=Q(i,j);
                M((i-1)*6+j,2)=RPM(i,j);
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
    if 1 %fim
        for i=1:4*length(Q)/17
            for j=1:6
                %M((i-1)*6+j,1)=Q(i+13*(length(Q)-1)/17,j);
                %M((i-1)*6+j,2)=RPM(i+13*(length(Q)-1)/17,j);
                %fprintf("%f,%f\n",Q(i,j),RPM(i,j));
            end
        end
    end
    if csv
        %strcsv=strcat('V4_200Hz_','G',num2str(gait),'-X',num2str(axf),'-Z',num2str(azf),'-Tr',num2str(ang_tr)','-P12','-Pc0.1','-T-1','.csv');
        %csvwrite(strcsv,M);
        csvwrite('Fim G0.3 Ax 28.937 Az 39.815 Ang_r 12.299',M);
    end    
end
j=length(Q); %ÚLTIMAS LINHAS DO CÓDIGO
for i=0:j-1 
    Time(i+1)=i*taxa;
end 