clear;
clc;


addpath(genpath('func'));
addpath(genpath('measure'));
addpath(genpath('update'));


load('datasets\ORL.mat');
Dataname = 'ORL';
Y = gt;
X{1}=X1';X{2}=X2';X{3}=X3';
num_view = length(X); 
fea = cell(num_view, 1);
for v = 1:num_view
    fea{v} = zscore(X{v}'); 
end

cls_num = length(unique(Y));  


lambda1_list=[1e-3];
lambda2_list=[1e-1];
d_list=[cls_num];%k,2k,3k
p_list=[0.4];


R = [];              
para = [];          
MaxAcc = 0;


for lambda1 = lambda1_list
    for lambda2 = lambda2_list
        for d = d_list
            for p = p_list
                fprintf('Start Running: %s, lambda1 = %f, lambda2 = %f,  d = %f, p = %f\n',Dataname, lambda1 , lambda2, d,p);
                tic;
             
                [Z] = DSRMVC(fea, lambda1, lambda2, d, p);
                time = toc;

                for iv = 1:10
                    % Perform spectral clustering
                    Predicted(iv,:) = SpectralClustering(Z, cls_num);

                    % Evaluate the clustering result
                    result(iv,:) = ClusteringMeasure(Y, Predicted(iv,:));

                end
                res_mean = mean(result);
                res_std = std(result);
          
                if res_mean(1) > MaxAcc
                    MaxAcc = res_mean(1);
                end
                fprintf('ACC = %.4f, NMI = %.4f, Purity = %.4f, Fscore = %.4f, Precision = %.4f, Recall = %.4f MaxACC = %.4f, time = %.2fs\n', ...
                    res_mean(1), res_mean(2), res_mean(3), res_mean(4),res_mean(5),res_mean(6),MaxAcc,  time);
            

            end
        end
    end
end






%

