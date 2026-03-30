function Ez = Get_Field_From_ABCoe(z_vec, ABCoe, d_layer, kz_layer)
    num_layer = length(kz_layer);    
    Ez = zeros(size(z_vec));
    
    % 1. 确定物理界面位置 (与你的 d_layer 长度对应)
    % z=0 是第一层和第二层的界面
    z_pts = zeros(num_layer + 1, 1);
    z_pts(1) = -1e6; % 第一层空气起始
    z_pts(2) = 0;    % 第二层起始
    for l = 2:num_layer-1
        z_pts(l+1) = z_pts(l) + d_layer(l-1);
    end
    z_pts(num_layer+1) = 1e6; % 最后一层结束
    
    % 2. 逐点计算场强
    for l = 1:num_layer
        mask = (z_vec >= z_pts(l)) & (z_vec < z_pts(l+1));
        if any(mask)
            A = ABCoe{l}(1);
            B = ABCoe{l}(2);
            kz = kz_layer(l);
            
            % --- 核心逻辑：确定参考点 z_ref ---
            if l < num_layer
                % 根据你的 M 矩阵定义，前 N-1 层系数定义在右界面
                z_ref = z_pts(l+1); 
            else
                % 最后一层没有 P 矩阵级联，系数定义在左界面
                z_ref = z_pts(l);
            end
            
            % 计算相对于参考点的位移 dz
            dz = z_vec(mask) - z_ref;
            
            % 计算标量场 (对于 TM, 这通常代表磁场 Hy 或电场平行分量)
            Ez(mask) = A * exp(1i * kz * dz) + B * exp(-1i * kz * dz);
        end
    end
end