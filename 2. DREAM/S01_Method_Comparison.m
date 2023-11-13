%% RENDOR improves the accuracy of Gene Regulatory Networks inference 

clear
clc
addpath(genpath(pwd));

%%
auc_before=zeros(1,9);
auc_after_ND=zeros(1,9);
auc_after_RENDOR=zeros(1,9);
auc_after_Silencer=zeros(1,9);
auc_after_ICM=zeros(1,9);

pr_before=zeros(1,9);
pr_after_ND=zeros(1,9);
pr_after_RENDOR=zeros(1,9);
pr_after_Silencer=zeros(1,9);
pr_after_ICM=zeros(1,9);

methods={'CLR';'Relavance';'ARACNE';'Pearson';'Spearman';'GENIE3';'TIGRESS';'Inferelator';'ANOVerence'};
m=4;

disp('******************* Evalution: *******************')

for j=1:9
    
    disp(['method name:',cell2mat(methods(j))])
    
    %% load data
    input_network=load_data(1,cell2mat(methods(j)));
    gold_file = './2. DREAM/DREAM Networks/DREAM5_NetworkInference_Edges_Network1.tsv';
    gold_edges = load_dream_network(gold_file);
    
    %% compute original performance
    prediction_before=change_network_format(input_network);
    [auc_before(j), pr_before(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_before);
    
    %% apply denoising algorithms
    input_network1 = (input_network-min(min(input_network)))./(max(max(input_network))-min(min(input_network)));
    output_network_ND = ND_regulatory(input_network1);
    output_network_RENDOR = RENDOR(input_network1,m);
    output_network_Silencer = Silencer(input_network);
    output_network_ICM = ICM(input_network);
    
    %% compute performance after RENDOR
    prediction_after_ND=change_network_format(output_network_ND);
    [auc_after_ND(j), pr_after_ND(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_after_ND);
    
    prediction_after_RENDOR=change_network_format(output_network_RENDOR);
    [auc_after_RENDOR(j), pr_after_RENDOR(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_after_RENDOR);
    
    prediction_after_Silencer=change_network_format(output_network_Silencer);
    [auc_after_Silencer(j), pr_after_Silencer(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_after_Silencer);
    
    prediction_after_ICM=change_network_format(output_network_ICM);
    [auc_after_ICM(j), pr_after_ICM(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_after_ICM);

    
end


auc_score=[auc_before; auc_after_RENDOR;auc_after_ND;auc_after_ICM;auc_after_Silencer];
pr_score=[pr_before;pr_after_RENDOR;pr_after_ND;pr_after_ICM;pr_after_Silencer];
writematrix(auc_score, './2. DREAM/Results/AUROC.csv')
writematrix(pr_score, './2. DREAM/Results/AUPR.csv')


