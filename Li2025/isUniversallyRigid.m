function flag = isUniversallyRigid(stress_mat, dim)
    if issymmetric(stress_mat)
        eigval = eig(stress_mat);
        isPSD = all(eigval >= 0);
        if isPSD
            if rank(stress_mat,1e-3) == length(stress_mat) - dim - 1
                flag = 1;
                disp("Framework is universally rigid.");
            else
                flag = 0;
                disp("Framework is not universally rigid.");
            end
        else
            flag = 0;
            disp("Framework is not PSD.");
        end
    else
        flag = 0;
        disp("Matrix is symmetric.");
    end

