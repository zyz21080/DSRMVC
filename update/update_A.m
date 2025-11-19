function A = update_A(B, P, C, E2, Y2, mu2,nV)


    [d, n] = size(B{1});    % 特征维度和样本数
    I = eye(n);             % 单位矩阵

    % 计算 M2 和 M2'
    M2 = I - C;             % (I - C)
    
    % 用 cell 存储每个视图的 A
    % A = cell(1, V);  % 初始化为 cell 数组

    temp1 = 0;
    temp2 = 0;
    for v = 1:nV
        % 定义左侧矩阵和右侧矩阵
        temp1 = temp1 + 2 * (B{v} * B{v}');
        A1 = temp1;
        temp2 = temp2 + 2*(B{v}*P{v}');

        % A1 = 2 * (B{v} * B{v}'); 
        B1 = mu2*I-mu2*C'-(mu2/2)*C-(mu2/2)*C'+mu2*C*C';
        C1 = temp2 + Y2*C'-Y2+ mu2*E2-mu2*E2*C';
        % 使用 Sylvester 方程求解 A
        
    end
    A = sylvester(A1, B1, C1 );
    A = max(A,0);
    % A = sylvester(temp1, B1, temp2 + mu2 * (E2-Y2/mu2)* M2' );  % 求解方程 AX + XA' = D
end