function E1 = update_E1(B, S, Y1, lambda1, mu1,nV)
% update_E1 求解 E1
%
% 输入：
%   B    : 1×V cell，每个 B{v} ∈ R^{d×n}
%   S    : 1×V cell，每个 S{v} ∈ R^{n×n}
%   Y1   : 1×V cell，每个 Y1{v} ∈ R^{d×n}
%   lambda1 : 标量，惩罚项参数
%   mu1 : 标量，增广拉格朗日参数
%
% 输出：
%   E1   : 1×V cell，每个 E1{v} ∈ R^{d×n}
  

    for v = 1:nV
        % 计算 E1 的解
        term =  mu1 * (B{v} - B{v} * S{v}) + Y1{v};
        E1{v} = term / (2 * lambda1 + mu1);  % 根据公式直接计算
    end
end
