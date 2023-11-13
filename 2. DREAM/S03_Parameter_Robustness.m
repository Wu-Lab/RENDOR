%% RENDOR improves the accuracy of Gene Regulatory Networks inference 

clear
clc
addpath(genpath(pwd));

%%
m = [1.2:0.1:2, 4:2:10];
methods={'CLR';'Relavance';'ARACNE';'Pearson';'Spearman';'GENIE3';'TIGRESS';'Inferelator';'ANOVerence'};
auc_before=zeros(1,length(methods));
auc_after_RENDOR=zeros(length(m),length(methods));
pr_before=zeros(1,length(methods));
pr_after_RENDOR=zeros(length(m),length(methods));


disp('******************* Evalution: *******************')

for j=1:length(methods)
    
    disp(['method name:',cell2mat(methods(j))])
    
    %% load data
    input_network=load_data(1,cell2mat(methods(j)));
    gold_file = '2. DREAM/DREAM Networks/DREAM5_NetworkInference_Edges_Network1.tsv';
    gold_edges = load_dream_network(gold_file);
    
    %% compute performance before RENDOR
    prediction_before=change_network_format(input_network);
    [auc_before(j), pr_before(j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_before);
    
    %% apply RENDOR algorithm 
    input_network1 = (input_network-min(min(input_network)))./(max(max(input_network))-min(min(input_network)));

    for i=1:length(m)
        output_network_RENDOR = RENDOR(input_network1,m(i));
        disp(m(i))
        %% compute performance after RENDOR
        prediction_after_RENDOR=change_network_format(output_network_RENDOR);
        [auc_after_RENDOR(i, j), pr_after_RENDOR(i, j)] = DREAM5_Challenge4_Evaluation(gold_edges, prediction_after_RENDOR);

    end
    
end


writematrix(auc_after_RENDOR, './2. DREAM/Results/para_auc_after_RENDOR.csv')
writematrix(pr_after_RENDOR, './2. DREAM/Results/para_pr_after_RENDOR.csv')
writematrix(auc_before, './2. DREAM/Results/para_auc_before.csv')
writematrix(pr_before, './2. DREAM/Results/para_pr_before.csv')
