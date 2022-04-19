classdef tensor
    methods(Static)
%% Hadamard Product

% This function computes the Hadarmard Product of two given matrices.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created: 2022

function C = mtx_prod_had(A,B)
    
    [ia,ja] = size(A);
    [ib,jb] = size(B);
    
    if (ia ~= ib) || (ja~=jb)
        disp('Invalid Matrices!')
        return;
    else
        C = A.*B;
        %C = zeros(ia,ja);
        %for i = 1:ia 
            %for j = 1:ja
                %C(i,j) = A(i,j)*B(i,j);
            %end
        %end
    end
end

%% Kronecker Product

% This function computes the Kronecker Product of two given matrices. 
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created:2022

function C = mtx_prod_kron(A,B)
    
    [ia,ja] = size(A);
    [ib,jb] = size(B);
    
    A = repelem(A,ib,jb);
    B = repmat(B,[ia ja]);
    C = A.*B;
    
end

%% Khatri-Rao Product

% This function computes the Khatri-Rao Product of two given matrices. 
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created:2022

function C = mtx_prod_kr(A,B)
    [ia,ja] = size(A);
    [ib,jb] = size(B);
    
    if (ja~=jb)
        disp('Invalid Matrices!')
        return;
    else
        C = zeros(ia*ib,ja);
        for j = 1:ja
            C(:,j) = tensor.mtx_prod_kron(A(:,j),B(:,j));
        end
    end
end

%% Least-Squares Khatri-Rao Factorization (LSKRF)

% This function computes the LSKRF of a given matrix.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created:2022

function [Ahat,Bhat] = LSKRF(C,ia,ib)
    [~, jc] = size(C);
    
    Ahat = complex(zeros(ia,jc),0);
    Bhat = complex(zeros(ib,jc),0);
    
    for j = 1:jc
        Cp = C(:,j);
        Cp = reshape(Cp, [ib ia]);
        [U,S,V] = svd(Cp);
        Ahat(:,j) = sqrt(S(1,1)).*conj(V(:,1));
        Bhat(:,j) = sqrt(S(1,1)).*U(:,1);
    end
end

%% Least-Square Kronecker Product Factorization (LSKronF)

% This function computes the LSKronF of a given matrix.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created: 2022

function [Ahat,Bhat] = LSKronF(C,ia,ja,ib,jb)
    [ic,jc] = size(C);
    
    I = (ic/ia) + zeros(1,ia);
    J = (jc/ja) + zeros(1,ja);
    blocks_of_C = mat2cell(C,I,J);
    
    k = 1;
    Chat = complex(zeros(ib*jb,ia*ja),0);
    for j = 1:ja
        for i = 1:ia
            vec_of_block = cell2mat(blocks_of_C(i,j));
            vec_of_block = vec_of_block(:);
            Chat(:,k) = vec_of_block;
            k = k + 1;
        end
    end
    
    [U,S,V] = svd(Chat);
    ahat = sqrt(S(1,1)).*conj(V(:,1));
    bhat = sqrt(S(1,1)).*U(:,1);
    Ahat = reshape(ahat,[ia ja]);
    Bhat = reshape(bhat, [ib jb]);
end
   
%% Unfolding

% This function computes the unfolding of a given tensor in its matrix.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created: 19/04/2022

function [A] = unfold(ten,mode)
    dim = size(ten);
    order = 1:numel(dim);
    order(mode) = [];
    order = [mode order];
    A = reshape(permute(ten,order), dim(mode), prod(dim)/dim(mode));
end

%% Folding

% This function computes the folding of a given matrix into its tensor.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created: 19/04/2022

% It's interesting to see how the dimmensions get swapped by the unfolding
% so the understanding of the code is clear.
function [ten] = fold(A,dim,mode)
    order = 1:numel(dim);
    order(mode) = [];
    order = [mode order];
    dim = dim(order);
    ten = reshape(A,dim);
    
    if mode == 1
        ten = permute(ten,order);
    else
        order = 1:numel(dim);
        for i = 2:mode
            order([i-1 i]) = order([i i-1]);
        end
        ten = permute(ten,order);
    end
end

%% N-mode Product

% This function computes n-mode product of a set of matrices and a tensor.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created: 20/04/2022

function [ten] = n_mod_prod(ten,matrices,modes)
    dim = size(ten);
    number = numel(matrices);
    if nargin < 3
       modes = 1:number; 
    end
    
    for i = modes
        ten = cell2mat(matrices(i))*tensor.unfold(ten,i);
        [aux,~] = size(cell2mat(matrices(i)));
        dim(i) = aux;
        ten = tensor.fold(ten,[dim],i);
    end
end

%% High Order Single Value Decomposition (HOSVD)

% This function computes the HOSVD of a given tensor.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created: 21/04/2022

function [S,U] = HOSVD(ten,ranks)
    number = numel(size(ten));
    for i = 1:number
       [aux,~,~] = svd(tensor.unfold(ten,i)); 
       aux = aux(1:ranks(i),:);
       U{i} = aux;
    end
    % Core tensor uses the hermitian operator.
    Ut = cellfun(@(x) conj(x),U,'UniformOutput',false); 
    S = tensor.n_mod_prod(ten,Ut);
    % The normal factors should be transposed.
    U = cellfun(@(x) x.',U,'UniformOutput',false);
end

%% High Order Orthogonal Iteration (HOOI)

% This function computes the HOOI of a given tensor.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created: 21/04/2022

function [S,U] = HOOI(ten)
    max_iter = 5;
    [~, U] = tensor.HOSVD(ten);
    number = numel(size(ten));
    for k = 1:max_iter
        for i = 1:number
            modes = 1:number;
            modes(i) = []; % It will skip this mode in the n_mod_prod.
            Un = tensor.n_mod_prod(ten,U,modes);
            [aux,~,~] = svd(tensor.unfold(Un,i));
            U{i} = aux';
        end
    end
    S = tensor.n_mod_prod(ten,U);
end

%% Alternate Least-Square (ALS) [NEEDS CORRECTION]

% This function computes the ALS of a given tensor.   
% Author: Kenneth B. dos A. Benicio <kenneth@gtel.ufc.br>
% Created: 2022

function [A,B,C,error] = ALS(ten,R)
    I = zeros(R,R,R); 
    for i = 1:R
        I(i,i,i) = 1; 
    end
    [ia,ib,ic] = size(ten);
    
    mode_1 = tensor.unfold(ten,1);
    mode_2 = tensor.unfold(ten,2);
    mode_3 = tensor.unfold(ten,3);

    A = randn(ia,R) + 1j*randn(ia,R);
    B = randn(ib,R) + 1j*randn(ib,R);
    C = randn(ic,R) + 1j*randn(ic,R);

    aux = 10000;
    error = zeros(1,aux);
    error(1) = ((norm((mode_1 - B*(tensor.mtx_prod_kr(A,C).')),'fro'))^2)/((norm(mode_1,'fro')^2));
    for i = 2:aux
        B = mode_1*pinv((tensor.mtx_prod_kr(A,C)).');
        C = mode_2*pinv((tensor.mtx_prod_kr(A,B)).');
        A = mode_3*pinv((tensor.mtx_prod_kr(C,B)).');
        error(i) = ((norm((mode_1 - B*(tensor.mtx_prod_kr(A,C).')),'fro'))^2)/((norm(mode_1,'fro')^2));
        if abs(error(i) - error(i-1)) < eps
            error = error(1:i);
            break;
        else
            continue;
        end
    end
end
    
    end
end