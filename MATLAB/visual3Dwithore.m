load('ABC0006_withore.mat')
Gamma1s_orin = Gamma1s_orin*200;
% Gamma1s_ore_orin=Gamma1s_ore_orin*2;
figure;
for(id = 1 : size(Gamma1s_orin,2))
hold on;
colorcur = [0 0 1];%rand(1,3);
plot3(Gamma1s_orin(1,id), Gamma1s_orin(2,id), Gamma1s_orin(3,id), ...
    'x','Color',colorcur,'LineWidth',1,'MarkerSize',3)
plot3([Gamma1s_orin(1,id)-Gamma1s_ore_orin(1,id), Gamma1s_orin(1,id)+Gamma1s_ore_orin(1,id)], ...
      [Gamma1s_orin(2,id)-Gamma1s_ore_orin(2,id), Gamma1s_orin(2,id)+Gamma1s_ore_orin(2,id)], ...
      [Gamma1s_orin(3,id)-Gamma1s_ore_orin(3,id), Gamma1s_orin(3,id)+Gamma1s_ore_orin(3,id)],...
    'Color',colorcur,'LineWidth',1);
end
% axis equal
% view(3)

load('ABC0006_withore1.mat')
Gamma1s_orin = Gamma1s_orin*200;
% Gamma1s_ore_orin=Gamma1s_ore_orin*2;
% figure;
for(id = 1 : size(Gamma1s_orin,2))
hold on;
colorcur = [0 1 0];%rand(1,3);
plot3(Gamma1s_orin(1,id), Gamma1s_orin(2,id), Gamma1s_orin(3,id), ...
    'x','Color',colorcur,'LineWidth',1,'MarkerSize',3)
plot3([Gamma1s_orin(1,id)-Gamma1s_ore_orin(1,id), Gamma1s_orin(1,id)+Gamma1s_ore_orin(1,id)], ...
      [Gamma1s_orin(2,id)-Gamma1s_ore_orin(2,id), Gamma1s_orin(2,id)+Gamma1s_ore_orin(2,id)], ...
      [Gamma1s_orin(3,id)-Gamma1s_ore_orin(3,id), Gamma1s_orin(3,id)+Gamma1s_ore_orin(3,id)],...
    'Color',colorcur,'LineWidth',1);
end
% axis equal
% view(3)

load('ABC0006_withore2.mat')
Gamma1s_orin = Gamma1s_orin*200;
% Gamma1s_ore_orin=Gamma1s_ore_orin*2;
% figure;
for(id = 1 : size(Gamma1s_orin,2))
hold on;
colorcur = [1 0 1];%rand(1,3);
plot3(Gamma1s_orin(1,id), Gamma1s_orin(2,id), Gamma1s_orin(3,id), ...
    'x','Color',colorcur,'LineWidth',1,'MarkerSize',3)
plot3([Gamma1s_orin(1,id)-Gamma1s_ore_orin(1,id), Gamma1s_orin(1,id)+Gamma1s_ore_orin(1,id)], ...
      [Gamma1s_orin(2,id)-Gamma1s_ore_orin(2,id), Gamma1s_orin(2,id)+Gamma1s_ore_orin(2,id)], ...
      [Gamma1s_orin(3,id)-Gamma1s_ore_orin(3,id), Gamma1s_orin(3,id)+Gamma1s_ore_orin(3,id)],...
    'Color',colorcur,'LineWidth',1);
end
axis equal
view(3)

load('ABC0006_4viewpairs.mat')
Gamma1s_orin_1 = Gamma1s_orin_1*200;
Gamma1s_ore_orin1_1=Gamma1s_ore_orin1_1*1;
figure;
for(id = 1 : size(Gamma1s_orin_1,2))
hold on;
colorcur = [0 0 1];%rand(1,3);
plot3(Gamma1s_orin_1(1,id), Gamma1s_orin_1(2,id), Gamma1s_orin_1(3,id), ...
    'x','Color',colorcur,'LineWidth',1,'MarkerSize',3)
plot3([Gamma1s_orin_1(1,id)-Gamma1s_ore_orin1_1(1,id), Gamma1s_orin_1(1,id)+Gamma1s_ore_orin1_1(1,id)], ...
      [Gamma1s_orin_1(2,id)-Gamma1s_ore_orin1_1(2,id), Gamma1s_orin_1(2,id)+Gamma1s_ore_orin1_1(2,id)], ...
      [Gamma1s_orin_1(3,id)-Gamma1s_ore_orin1_1(3,id), Gamma1s_orin_1(3,id)+Gamma1s_ore_orin1_1(3,id)],...
    'Color',colorcur,'LineWidth',1);
end
% axis equal
% view(3)

Gamma1s_orin_2 = Gamma1s_orin_2*200;
Gamma1s_ore_orin1_2=Gamma1s_ore_orin1_2*1;
% figure;
for(id = 1 : size(Gamma1s_orin_2,2))
hold on;
colorcur = [0 1 0];%rand(1,3);
plot3(Gamma1s_orin_2(1,id), Gamma1s_orin_2(2,id), Gamma1s_orin_2(3,id), ...
    'x','Color',colorcur,'LineWidth',1,'MarkerSize',3)
plot3([Gamma1s_orin_2(1,id)-Gamma1s_ore_orin1_2(1,id), Gamma1s_orin_2(1,id)+Gamma1s_ore_orin1_2(1,id)], ...
      [Gamma1s_orin_2(2,id)-Gamma1s_ore_orin1_2(2,id), Gamma1s_orin_2(2,id)+Gamma1s_ore_orin1_2(2,id)], ...
      [Gamma1s_orin_2(3,id)-Gamma1s_ore_orin1_2(3,id), Gamma1s_orin_2(3,id)+Gamma1s_ore_orin1_2(3,id)],...
    'Color',colorcur,'LineWidth',1);
end
% axis equal
% view(3)

Gamma1s_orin_3 = Gamma1s_orin_3*200;
Gamma1s_ore_orin1_3=Gamma1s_ore_orin1_3*1;
% figure;
for(id = 1 : size(Gamma1s_orin_3,2))
hold on;
colorcur = [1 0 1];%rand(1,3);
plot3(Gamma1s_orin_3(1,id), Gamma1s_orin_3(2,id), Gamma1s_orin_3(3,id), ...
    'x','Color',colorcur,'LineWidth',1,'MarkerSize',3)
plot3([Gamma1s_orin_3(1,id)-Gamma1s_ore_orin1_3(1,id), Gamma1s_orin_3(1,id)+Gamma1s_ore_orin1_3(1,id)], ...
      [Gamma1s_orin_3(2,id)-Gamma1s_ore_orin1_3(2,id), Gamma1s_orin_3(2,id)+Gamma1s_ore_orin1_3(2,id)], ...
      [Gamma1s_orin_3(3,id)-Gamma1s_ore_orin1_3(3,id), Gamma1s_orin_3(3,id)+Gamma1s_ore_orin1_3(3,id)],...
    'Color',colorcur,'LineWidth',1);
end
% axis equal
% view(3)

Gamma1s_orin_4 = Gamma1s_orin_4*200;
Gamma1s_ore_orin1_4=Gamma1s_ore_orin1_4*1;
% figure;
for(id = 1 : size(Gamma1s_orin_4,2))
hold on;
colorcur = [0.9290 0.6940 0.1250];%rand(1,3);
plot3(Gamma1s_orin_4(1,id), Gamma1s_orin_4(2,id), Gamma1s_orin_4(3,id), ...
    'x','Color',colorcur,'LineWidth',1,'MarkerSize',3)
plot3([Gamma1s_orin_4(1,id)-Gamma1s_ore_orin1_4(1,id), Gamma1s_orin_4(1,id)+Gamma1s_ore_orin1_4(1,id)], ...
      [Gamma1s_orin_4(2,id)-Gamma1s_ore_orin1_4(2,id), Gamma1s_orin_4(2,id)+Gamma1s_ore_orin1_4(2,id)], ...
      [Gamma1s_orin_4(3,id)-Gamma1s_ore_orin1_4(3,id), Gamma1s_orin_4(3,id)+Gamma1s_ore_orin1_4(3,id)],...
    'Color',colorcur,'LineWidth',1);
end
axis equal
view(3)