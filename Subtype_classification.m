%% landscape 绘制Accuracy, Precision, Recall,f1-measure,roc曲线及计算AUC
% 导入数据为F:\SPA pythonFile\TCGA_LUAD_landscape\subtype\LUAD_entorpy_subtype.xls
load('subtype.mat');
min_entropy = min(entropy);
max_entropy = max(entropy);

point_num = 0;
for threshold = min_entropy:max_entropy+1 
    point_num = point_num + 1;
end

x = zeros(point_num, 1)   
TPR = zeros(point_num,1);
FPR = zeros(point_num,1);
Acc = zeros(point_num,1);
Precision = zeros(point_num,1);
Recall = zeros(point_num,1);
f_measure = zeros(point_num,1);
i = 1;
for threshold = min_entropy:max_entropy+1 
    TP = 0;
    FP = 0;
    FN = 0;
    TN = 0;
    x(i) = threshold;
    prodict_label = zeros(length(entropy),1);
    for j = 1:length(entropy)
        if entropy(j) >= threshold
            prodict_label(j) = 1;
        else
            prodict_label(j) = -1;
        end
    end
    for k = 1:length(entropy)
        if true_label(k) == 1 && prodict_label(k) == 1
            TP = TP + 1;
        elseif true_label(k) == 1 && prodict_label(k) == -1
            FN = FN + 1;
        elseif true_label(k) == -1 && prodict_label(k) == 1
            FP = FP + 1;
        else
            TN = TN + 1;
        end
    end
    TPR(i) = TP / (TP + FN);
    FPR(i) = FP / (FP + TN);
    Acc(i) = (TP + TN) / (TP + TN + FP +FN);
    Recall(i) = TP / (TP + FN);
    Precision(i) = TP / (TP + FP);
    f_measure(i) = 2*TP / (2*TP + FP + FN);
    i = i + 1;
end

AUC = 0;
for i = 1 : point_num - 1
    AUC = AUC + (FPR(i) - FPR(i + 1)) * TPR(i);
end

figure
plot(x, Acc)
hold on 
plot(x, Precision)
hold on
plot(x, Recall)
hold on
plot(x, f_measure)
hold on
xlabel('Threshold Value')
legend_names = {'Accuracy', 'Precision','Recall','f1-measure'};
legend(legend_names,'Location','northwest')

figure
plot(FPR,TPR);
hold on
x=0:0.1:1;
y=x;
plot(x,y)
title('ROC Curves')
xlabel('False Positive Rate')
ylabel('True Positive Rate')
legend(legend_names,'Location','southeast')

