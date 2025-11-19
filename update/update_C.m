function C = update_C(A, E2, Y2, Q, mu2, rho, J,nV)
   
    [d, n] = size(A);    % 特征维度和样本数
    I = eye(n);             % 单位矩阵

    % 计算 A' A
    AtA = A' * A;  % r×r

    % 计算第一项：(mu2/rho * A' A + I)
    left_term = (mu2 / rho) * AtA + eye(n);

    % 计算第二项：(mu2/rho * A' A - mu2 * A' * E2 + A' * Y2 - Q)
    right_term = (mu2 * AtA - mu2 * A' * E2 + A' * Y2 - Q{nV+1})/rho;

    % 计算 C
    C = inv(left_term)*(right_term + J{nV+1}) ;  % 求解并加上 J
end
