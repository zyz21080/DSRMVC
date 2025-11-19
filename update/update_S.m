function S = update_S(B, E1, Y1, Q, mu1, rho, J,nV)
% update_S 求解 S^{(v)}
%
% 输入：
%   B    : 1×V cell，每个 B{v} ∈ R^{d×n}
%   E1   : 1×V cell，每个 E1{v} ∈ R^{d×n}
%   Y1   : 1×V cell，每个 Y1{v} ∈ R^{d×n}
%   Q    : 1×V cell，每个 Q{v} ∈ R^{d×n}
%   mu1  : 标量，增广拉格朗日参数
%   rho  : 标量，控制增广项
%   J    : 1×V cell，每个 J{v} ∈ R^{d×n}


    [d, n] = size(B{1});    % 特征维度和样本数
    I = eye(n);             % 单位矩阵

    for v = 1:nV
        
        % 计算左侧矩阵
        left_term = (mu1 / rho) * (B{v}' * B{v}) + I;

        % 计算右侧矩阵
        right_term = (mu1* (B{v}' * B{v}) - mu1 * (B{v}' * E1{v}) + (B{v}' * Y1{v}) - Q{v})/rho;

        % 计算 S^{(v)}
        S{v} = inv(left_term)*(right_term + J{v});
    end
end
