function [Z] = DSRMVC(X,lambda1,lambda2,d,p)
    nV = length(X);
    [n] = size(X{1},2);
    Isconverg = 0;epson = 1e-5;
    iter = 0;
    maxiter=200;
    converge_Z1=[];
    converge_Z2=[];
    converge_Z_J=[];


    mu1       = 1e-3;
    mu2       = 1e-3;
    rho       = 1e-3;
    max_mu1  =  1e10;
    max_mu2   = 1e10;
    max_rho   = 1e10;
    sX = [n, n, nV+1];
    A  = zeros(d, n);
    C  = zeros(n, n);
    weight_vector = ones(1,nV+1)';
    for i=1:nV
        X{i} = X{i}./repmat(sqrt(sum(X{i}.^2,1)),size(X{i},1),1);
        % T{i} = constructW_PKN(X{i},40);  
        P{i} = SPPMI(constructW_PKN(X{i},40), 3); 
    end

   
    for v = 1:nV
        B{v}   = eye(d, n);      % B^{(v)} ≥ 0
        S{v}   = zeros(n, n);
        E1{v}  = zeros(d, n);
        Y1{v}  = zeros(d, n);
    end
    for v = 1:nV+1
        J{v} = zeros(n,n);
        Q{v} = zeros(n,n);
    end

    % —— 全局噪声项 & 乘子 ——
    E2 = zeros(d, n);
    Y2 = zeros(d, n);
   
    while(Isconverg == 0)
   
        % === Subproblem A ===
        A = update_A(B, P, C, E2, Y2, mu2,nV);

        % === Subproblem B ===
        B = update_B(A, P, S, E1, Y1, mu1,nV);

        % === Subproblem S ===
        S = update_S(B, E1, Y1, Q, mu1, rho, J,nV);

        % === Subproblem C ===
        C = update_C(A, E2, Y2, Q, mu2, rho, J,nV);

        % === Subproblem E1 ===
        E1 = update_E1(B, S, Y1, lambda1, mu1,nV);

        % === Subproblem E2 ===
        E2 = update_E2(A, C, Y2, lambda1, mu2,nV);
        % === Subproblem tensorJ ===
        Z_tensor = cat(3, S{:,:},C);
        Q_tensor = cat(3, Q{:,:});
        z = Z_tensor(:);
        q = Q_tensor(:);
        % [j,objv] = wshrinkObj(z + 1/rho*q,lambda2./rho,sX, 0,3);
        [j] = wshrinkObj_weight_lp(z + 1/rho*q, (weight_vector*lambda2)./rho,sX, 0,3,p);
        J_tensor = reshape(j, sX);
       
        for iv = 1 : nV
         S{iv} = J_tensor(:,:,iv);
        end
         C= J_tensor(:,:,nV+1);
       
         % === Update Lagrange multiplier ===
        for v=1:nV
            Y1{v}=Y1{v} + mu1 * (B{v} - B{v} * S{v} - E1{v});
        end
        Y2 = Y2 + mu2*(A - A*C - E2);
        
        for iv = 1 : nV+1
        Q_tensor(:,:,iv) = Q_tensor(:,:,iv) + rho*(Z_tensor(:,:,iv)-J_tensor(:,:,iv));
        Q{iv} = Q_tensor(:,:,iv);
        end
    
        % === Update mu ===
        mu1 = min(mu1*2, max_mu1);
        mu2 = min(mu2*2, max_mu2);
        rho = min(rho*2, max_rho);
        % === Convergence check ===
        max_Z1=0;
        max_Z2=0;
        max_Z_J=0;
        Isconverg = 1;
        for k=1:nV
            if (norm(B{k}-B{k}*S{k}-E1{k},inf)>epson)
                history.norm_Z1 = norm(B{k}-B{k}*S{k}-E1{k},inf);
                % fprintf("history.norm_Z1=%f\n,",history.norm_Z1);
                Isconverg = 0;
                max_Z1=max(max_Z1,history.norm_Z1);
            end

            J{k+1} = J_tensor(:,:,k+1);
            Z{k+1} = Z_tensor(:,:,k+1);
            Q_tensor = reshape(q, sX);
            Q{k+1} = Q_tensor(:,:,k+1);
            if (norm(Z{k+1}-J{k+1},inf)>epson)
            history.norm_Z_J = norm(Z{k+1}-J{k+1},inf);
            % fprintf("history.norm_Z_J=%f\n,",history.norm_Z_J);
            Isconverg = 0;
            max_Z_J=max(max_Z_J, history.norm_Z_J);
            end
        end
        if (norm(A-A*C-E2,inf)>epson)
            history.norm_Z2 = norm(A-A*C-E2,inf);
            % fprintf("history.norm_Z2=%f\n,",history.norm_Z2);
            Isconverg = 0;
            max_Z2=max(max_Z2,history.norm_Z2);
        end
    
        if (iter>maxiter)
            Isconverg  = 1;
        end
        iter = iter + 1;
    end
TMP = 0;
for v=1:nV
    TMP = TMP + abs(S{v}) + abs(S{v}');
end
Z = (TMP/nV + abs(C) + abs(C'))/2;

fprintf("迭代了%d次后收敛\n",iter);

