% P_edge = [];
% R_edge = [];

lineColor = [0,1,0; ...
             0,0,1; ...
             0,1,1; ...
             1,0,1; ...
             0.8500,0.3250,0.0980; ...
             0.9290,0.6940,0.1250; ...
             0.4940,0.1840,0.5560; ...
             0.4660,0.6740,0.1880; ...
             0.3010,0.7450,0.9330; ...
             0.6350,0.0780,0.1840];

figure
for(i = 1:1: size(P_edge,1))
    hold on
%     endidxall = find(R_edge(i,:)~=0);
% p = plot(R_edge(i,1:endidxall(end)), P_edge(i,1:endidxall(end)),'LineWidth',1,'Color',lineColor(i,:));
% p.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Presicion',P_edge(i,1:endidxall(end)));
% p.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Recall',R_edge(i,1:endidxall(end)));
stidxall = find(R_edge(i,:)~=0);
p = plot(R_edge(i,stidxall), P_edge(i,stidxall),'LineWidth',1,'Color',lineColor(i,:),'Marker','square');
p.DataTipTemplate.DataTipRows(1) = dataTipTextRow('Presicion',P_edge(i,stidxall));
p.DataTipTemplate.DataTipRows(2) = dataTipTextRow('Recall',R_edge(i,stidxall));
% p.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Δθ',[1,2,5,10,15,30]);
% p.DataTipTemplate.DataTipRows(3) = dataTipTextRow('N support',[4,6,8,10,12,16]);
% p.DataTipTemplate.DataTipRows(3) = dataTipTextRow('Δ',[0.05,0.08, 0.1,0.3,0.6,1.2]);
end
xlabel ('Recall','FontSize',15)
ylabel ('Precision','FontSize',15)
% title  ('Precision-Recall Curve (changing Δθ)','FontSize',15)
% title  ('Precision-Recall Curve (changing N support)','FontSize',15)
% title  ('Precision-Recall Curve (changing Δ)','FontSize',15)
% legend ({'Δ=0.05','Δ=0.08','Δ=0.1','Δ=0.3','Δ=0.6','Δ=1.2'},'FontSize',10)
legend ({'Δθ=1°','Δθ=2°','Δθ=5°','Δθ=10°','Δθ=15°','Δθ=30°'},'FontSize',10)
% legend ({'N=4','N=6','N=8','N=10','N=12','N=14'},'FontSize',10)
xlim([0,1])
ylim([0,1.01])