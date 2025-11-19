function B = update_B(A, P, S, E1, Y1, mu1,nV)


    [n, n] = size(S{1});  % 数据的维度

    % 用 cell 存储每个视图的 B
    % B = cell(1, length(S));

    % 循环每个视图
    for v = 1:nV
        I = eye(n);
        A1 = 2*(A*A');
        B1= mu1*I-mu1*S{v}'-(mu1/2)*S{v}-(mu1/2)*S{v}'+mu1*S{v}*S{v}';
        C1= 2 * A * P{v} + Y1{v}*S{v}- Y1{v} + mu1*E1{v}-mu1*E1{v}*S{v}';
        B{v} = sylvester(A1, B1, C1); 
        B{v} = max(B{v},0);
    end
end
