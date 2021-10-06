%Mutation
P_A=Pp;
F_V_A=F_V;
% k diferente de pmax significa que o melhor indivduo no sofre
% mutação
%k =1 todos sofrem mutaço
for k=1:length(Pp) 
    if 1 %k ~=pmax 
        Pp(k,1)=Pp(k,1)-0.05+0.1*rand(1);  %angulo do torso
        Pp(k,2)=Pp(k,2)-0.5+rand(1);%angulo de inclinação
        Pp(k,3)=Pp(k,3)-0.5+rand(1);  %angulo do torso
        Pp(k,4)=Pp(k,4)-5+10*rand(1); %amplitude horizontal
        Pp(k,5)=Pp(k,5)-5+10*rand(1); %amplitude vertical
    end
end