function E2 = update_E2(A, C, Y2, lambda1, mu2,nV)
% update_E2 求解 E2
%
% 输入：
%   A    : d×r 矩阵，公共子空间表示
%   C    : r×r 矩阵
%   Y2   : d×r 矩阵，拉格朗日乘子
%   lambda1 : 标量，惩罚项参数
%   mu2 : 标量，增广拉格朗日参数
%
% 输出：
%   E2   : d×r 矩阵，误差项 E2

    % 计算 E2 的解
    E2 =  (mu2 * A - mu2 * A * C + Y2)/(2 * lambda1 + mu2) ;
end
