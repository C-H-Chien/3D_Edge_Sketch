load fisheriris
pred = meas(51:end,1:2);
resp = (1:100)'>50;  % Versicolor = 0, virginica = 1
mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
scores = mdl.Fitted.Probability;
[X,Y,T,AUC] = perfcurve(species(51:end,:),scores,'virginica', 'XCrit', 'tpr', 'YCrit', 'prec');
AUC
plot(X, Y);
xlabel('Recall')
ylabel('Precision')
xlim([0, 1])
ylim([0, 1])

repsondt = 0.002;
