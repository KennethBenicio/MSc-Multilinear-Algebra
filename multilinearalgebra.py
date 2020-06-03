# -------------------------------------------------
# Module containing the implemented functions  
# during the Tensor Algebra course. 
# These functions are composing my own tensor 
# toolbox.
# -------------------------------------------------
## Author: Kenneth Benício
## Version: 0.0.0
## Email: Kennethbrenner4242@gmail.com
## Status: in development

# Imports
import numpy as np


# Usefull functions

def vec(matrix):
    return matrix.flatten(order='F').reshape(-1, 1)


def unvec(vector, nrow, ncol):
    return vector.reshape(nrow, ncol, order='F')


def extract_block(matrix, shape_A, shape_B):
    shape_A = np.array(shape_A)
    shape_B = np.array(shape_B)
    ratio = int(matrix.shape[1] / shape_B[1])
    chunks = np.split(matrix, ratio, axis=1)

    tensor = np.array(chunks).reshape(shape_A.prod(), shape_B[0], shape_B[1])

    matrix_approx = np.zeros([shape_B.prod(), shape_A.prod()], dtype='complex')
    for s in range(shape_A.prod()):
        matrix_approx[:, s] = tensor[s].flatten(order='F')

    return matrix_approx


def extract_block_mlskronf(matrix, shapes):
    nrow_a, ncol_a = shapes[0]
    nrow_b, ncol_b = [int(x / y) for x, y in zip(matrix.shape, shapes[0])]

    split = np.array(np.hsplit(matrix, ncol_a))

    aux = split.reshape(nrow_a * ncol_a, nrow_b, ncol_b)

    nrow_a, ncol_a = shapes[1]
    nrow_b, ncol_b = [int(x / y) for x, y in zip((nrow_b, ncol_b), shapes[1])]

    if np.iscomplexobj(matrix):
        out_x = np.zeros((nrow_b * ncol_b * nrow_a * ncol_a, np.prod(shapes[0])),
                         dtype=np.complex_)
    else:
        out_x = np.zeros((nrow_b * ncol_b * nrow_a * ncol_a, np.prod(shapes[0])))

    for s in range(aux.shape[0]):
        split_aux = np.array(np.hsplit(aux[s], ncol_a))
        split_aux = split_aux.reshape(nrow_a * ncol_a, nrow_b, ncol_b)
        out_x[:, [s]] = vec(split_aux.T.reshape(nrow_b * ncol_b, nrow_a * ncol_a))

    return out_x


def mode(matrix, tensor):
    matrix = matrix.shape
    tensor = tensor.shape

    for i in range(len(matrix)):
        for j in range(len(tensor)):
            if matrix[i] == tensor[j]:
                return j


def normalized_mean_square_error(tensor, tensor_hat):
    nmse = (np.linalg.norm(tensor - tensor_hat)) ** 2 / (np.linalg.norm(tensor)) ** 2

    return nmse


def tensor_norm(tensor):
    return np.sqrt(np.sum(tensor ** 2))


def build_cpd(shapes, rank):
    A = np.random.randn(shapes[1], rank)
    B = np.random.randn(shapes[2], rank)
    C = np.random.randn(shapes[0], rank)

    cpd_tensor = np.zeros(shapes)
    for i in range(0, rank):
        a = (A[:, i])[:, None]
        b = (B[:, i])[:, None]
        c = (C[:, i])[:, None]

        tensor = (a @ b.T) * (c[:, :, None])
        cpd_tensor = cpd_tensor + tensor

    return cpd_tensor


# Matrix products

def hadamard(list_of_arrays):
    matrix = np.ones(list_of_arrays[0].shape)
    Size = len(list_of_arrays)
    for i in range(0, Size):
        matrix = matrix * list_of_arrays[i]

    return matrix


def kron(list_of_arrays):
    shapes = []
    for i in range(0, len(list_of_arrays)):
        shape_of_array = list_of_arrays[i].shape
        shapes.append(shape_of_array)

    matrix = list_of_arrays[0]
    for i in range(1, len(list_of_arrays)):
        matrix = matrix[:, np.newaxis, :, np.newaxis] * list_of_arrays[i][np.newaxis, :, np.newaxis, :]

        shapes_aux = np.ones([np.size(shapes[0])], dtype='int8')
        for j in range(i - 1, i + 1):
            for k in range(0, np.size(shapes[0])):
                shapes_aux[k] = shapes_aux[k] * shapes[j][k]

        matrix = matrix.reshape(shapes_aux)
        shapes[i] = matrix.shape

    return matrix


def khatri(list_of_arrays):
    shapes = []
    for i in range(0, len(list_of_arrays)):
        shape_of_array = list_of_arrays[i].shape
        shapes.append(shape_of_array)

    shapes_aux = np.ones([np.size(shapes[0])], dtype='int8')
    for j in range(0, len(list_of_arrays)):
        for k in range(0, np.size(shapes[0])):
            shapes_aux[k] = shapes_aux[k] * shapes[j][k]

    matrix = np.ones([shapes_aux[0], list_of_arrays[0].shape[1]])

    for j in range(0, list_of_arrays[0].shape[1]):
        list_to_kronecker = []
        for i in range(0, len(list_of_arrays)):
            list_to_kronecker.append((list_of_arrays[i][:, j])[:, None])

    for k in range(0, list_of_arrays[0].shape[1]):
        matrix[:, k] = (kron(list_to_kronecker)).flatten()

    return matrix


# Special Matrix Factorizations

def LSKRF(matrix, I, J, R):
    # X: A matrix created by the relation X = A ⋄ B ∈ C^{IJ×R}, where A ∈ C^{I×R} and B ∈ C^{J×R}.
    # (I,J,R): Dimensions of A and B as mentioned above.

    # Creating a matrix to alocate the values of stimated X.
    matrix_hat = np.zeros([I * J, R]) + 1j * np.zeros([I * J, R])

    for i in range(0, R):
        # Each column of X will be rearranged into a matrix Xp ∈ ℂ^{J×I}.
        matrix_p = (matrix[:, i]).reshape(J, I, order='F')

        # In this line the SVD of the new matrix will be calculated.
        [U_p, S_p, V_p] = np.linalg.svd(matrix_p)

        # Now the stimations of the a_hat and b_hat i-th columns will be made.
        a_hat = np.sqrt(S_p[0]) * ((V_p[0, :]))
        b_hat = np.sqrt(S_p[0]) * U_p[:, 0]

        # Finally, the kronecker between the stimations of the i-th columns of a_hat and b_hat
        # will be calculated and alocated to the respective column of the stimation of X.
        matrix_hat[:, i] = (np.outer(b_hat, a_hat)).reshape(I * J, order='F')

    # This line will calculate the Normalized Mean Square Error between the stimation X_hat and the know signal X.
    nmse = normalized_mean_square_error(matrix, matrix_hat)

    return matrix_hat, nmse


def LSKronF(matrix, I, P, J, Q):
    # X: The uncorrupted signal X = A ⊗ B ∈ C^{IJ×PQ}, where A ∈ C^{I×P} and B ∈ C^{J×Q}.
    # (I,P,J,Q): Dimensions of A and B as mentioned above.

    # This line will call a function to construct the new matrix X
    # as X_approx = vec(B)vec(A)^T that will be used in the next computations.
    matrix_approx = extract_chunk(matrix, [I, P], [J, Q])
    # Obs:Is worth note that, except for the dimensions of A and B, the previously knowledge of the data
    # in A and B are obviously unknow.

    # In this line the SVD of the new matrix will be calculated.
    [U_p, S_p, V_p] = np.linalg.svd(matrix_approx)

    # Now the stimations of vec(A) and vec(B) will be made.
    A_hat = (np.sqrt(S_p[0]) * ((V_p[0, :])))[:, None]
    B_hat = (np.sqrt(S_p[0]) * U_p[:, 0])[:, None]

    # And here the unvec operation is made to obtain the final stimations of A and B.
    A_hat = A_hat.reshape(I, P, order='F')
    B_hat = B_hat.reshape(J, Q, order='F')

    matrix_hat = tl.tenalg.kronecker([A_hat, B_hat])

    # This line will calculate the Normalized Mean Square Error between the stimation X_hat and the know signal X.
    nmse = normalized_mean_square_error(matrix, matrix_hat)

    return matrix_hat, nmse


# Tensor Operations

def unfold(tensor, n):
    tensor = np.moveaxis(tensor, n, 0)
    tensor_unfolding = tensor.reshape(tensor.shape[0], -1)

    return tensor_unfolding


def fold(tensor_unfolding, tensor_shape, n):
    # Transforming the shape of tensor tuple into a list for easy manipulation.
    shape = list(tensor_shape)
    # Extracting the external dimension that is presented in the unfolding tensor as the number of rows.
    n_dimension = shape.pop(n)
    # Inserting the previously dimension at the begining of the shape vector so this way we have a dinamic reshape
    # that will change in accord with the unfolding mode.
    shape.insert(0, n_dimension)

    # Reorganizing the unfolded tensor as a tensor.
    tensor = tensor_unfolding.reshape(shape)

    # Moving back the axis that was changed at the unfolding function.
    tensor = np.moveaxis(tensor, 0, n)

    return tensor


# DISCLAIMER: I'm not the author of _kolda_reorder, kolda_unfold and kolda_fold.

def _kolda_reorder(ndim, mode):
    indices = list(range(ndim))
    element = indices.pop(mode)

    return ([element] + indices[::-1])


def kolda_unfold(tensor, mode):
    matrix = np.transpose(tensor, _kolda_reorder(tensor.ndim, mode)).reshape((tensor.shape[mode], -1))

    return matrix


def kolda_fold(matrix, mode, shape):
    unfolded_indices = _kolda_reorder(len(shape), mode)
    original_shape = [shape[i] for i in unfolded_indices]
    matrix = matrix.reshape(original_shape)

    folded_indices = list(range(len(shape) - 1, 0, -1))
    folded_indices.insert(mode, 0)
    tensor = np.transpose(matrix, folded_indices)
    return tensor


def ten_mat_prod(tensor, matrix, mode):
    shape = list(tensor.shape)
    shape.pop(mode)
    shape.insert(mode, matrix.shape[0])

    tensor = matrix @ unfold(tensor, mode)
    tensor = fold(tensor, shape, mode)

    return tensor


def ten_mat_multiprod(tensor, list_of_matrices, *modes):
    if len(modes) == 0:
        modes = np.arange(0, list_of_matrices.shape[0])
        k = np.array(modes)[-2:]
        modes = np.append(k, (np.array(modes)[0:-2:]))

    shape = list(tensor.shape)
    for i in range(0, list_of_matrices.shape[0]):
        Z = list_of_matrices[i]
        tensor = Z @ (unfold(tensor, modes[i]))
        shape.pop(modes[i])
        shape.insert(modes[i], Z.shape[0])
        tensor = fold(tensor, shape, modes[i])
        shape = list(tensor.shape)

    return tensor


def ten_mat_vec(tensor, list_of_vectors, *modes):
    if len(modes) == 0:
        modes = np.arange(0, list_of_matrices.shape[0])
        k = np.array(modes)[-2:]
        modes = np.append(k, (np.array(modes)[0:-2:]))

    shape = list(tensor.shape)
    for i in range(0, list_of_vectors.shape[0]):
        Z = list_of_vectors[i]
        tensor = Z @ (unfold(tensor, modes[i]))
        shape.pop(modes[i])
        shape.insert(modes[i], Z.shape[0])
        tensor = fold(tensor, shape, modes[i])
        shape = list(tensor.shape)

    return tensor


# Tensor decompositions

def HOSVD(tensor, *ranks):
    # Full-rank HOSVD
    if len(ranks) == 0:
        U = []
        for i in range(0, tensor.ndim):
            [u, _, _] = np.linalg.svd(unfold(tensor, i))
            u = u.conj().T
            U.append(u)

        if len(np.array(U).shape) == 1:
            k = np.array(U)[-2:]
            U = np.append(k, (np.array(U)[0:-2:]))

        else:
            k = np.array(U)[-2:]
            U = np.append(k, (np.array(U)[0:-2:])).reshape(np.array(U).shape)

        S = ten_mat_multiprod(tensor, np.array(U))

        for i in range(0, tensor.ndim):
            U[i] = U[i].conj().T

    # Truncated HOSVD
    else:
        U = []
        for i in range(0, tensor.ndim):
            [u, s, _] = np.linalg.svd(unfold(tensor, i))
            u = u[:, 0:ranks[i]]
            u = u.conj().T
            U.append(u)

        if len(np.array(U).shape) == 1:
            k = np.array(U)[-2:]
            U = np.append(k, (np.array(U)[0:-2:]))

        else:
            k = np.array(U)[-2:]
            U = np.append(k, (np.array(U)[0:-2:])).reshape(np.array(U).shape)

        S = ten_mat_multiprod(tensor, np.array(U))

        U_aux = []
        for i in range(0, tensor.ndim):
            u = U[i].conj().T
            U_aux.append(u)

        U = U_aux

    return S, U, S.shape


def HOSVD_epsilon(tensor, epsilon):
    U = []
    for i in range(0, tensor.ndim):
        [u, s, _] = np.linalg.svd(unfold(tensor, i))
        s[s < epsilon] = 0
        s = [i for i in s if i != 0]
        s = len(s)
        u = u[:, 0:s]
        u = u.conj().T
        U.append(u)

    if len(np.array(U).shape) == 1:
        k = np.array(U)[-2:]
        U = np.append(k, (np.array(U)[0:-2:]))

    else:
        k = np.array(U)[-2:]
        U = np.append(k, (np.array(U)[0:-2:])).reshape(np.array(U).shape)

    S = ten_mat_multiprod(tensor, np.array(U))

    for i in range(0, tensor.ndim):
        U[i] = U[i].conj().T

    return S, U, S.shape


def HOOI(tensor, epsilon, *ranks):
    # Full-rank HOOI
    if len(ranks) == 0:
        [S_init, U_init, _] = HOSVD(tensor)

        for k in range(0, 10):
            U = []

            for i in range(0, tensor.ndim):
                # Matrix Selection
                aux = np.ones(tensor.ndim, dtype=bool)
                aux[i] = False
                U_aux = np.asarray(U_init)
                U_aux = U_aux[aux]

                # Creating the list of U matrices
                modes = np.zeros([tensor.ndim - 1], dtype='int8')
                for j in range(0, tensor.ndim - 1):
                    modes[j] = multilinearalgebra.mode(U_aux[j], tensor)

                u = multilinearalgebra.ten_mat_multiprod(tensor, np.array(U_aux), *modes)
                [u, _, _] = np.linalg.svd(multilinearalgebra.unfold(u, i))
                u = u.conj().T
                U.append(u)

            # Creating the core tensor
            S = multilinearalgebra.ten_mat_multiprod(tensor, np.array(U), *np.arange(0, tensor.ndim))

            # Convergence
            if multilinearalgebra.normalized_mean_square_error(S, S_init) > epsilon:
                # print('NMSE Error for the iteration',k + 1,'.')
                # print(multilinearalgebra.normalized_mean_square_error(S,S_init))
                S_init = S
                U_init = U

            else:
                for i in range(0, tensor.ndim):
                    U[i] = U[i].conj().T

                break

    # Truncated HOOI
    else:
        [S_init, U_init, _] = HOSVD(tensor, *ranks)

        for k in range(0, 10):
            U = []

            for i in range(0, tensor.ndim):
                # Matrix Selection
                aux = np.ones(tensor.ndim, dtype=bool)
                aux[i] = False
                U_aux = np.asarray(U_init)
                U_aux = U_aux[aux]

                # Creating the list of U matrices
                modes = np.zeros([tensor.ndim - 1], dtype='int8')
                for j in range(0, tensor.ndim - 1):
                    modes[j] = multilinearalgebra.mode(U_aux[j], tensor)

                u = multilinearalgebra.ten_mat_multiprod(tensor, np.array(U_aux), *modes)
                [u, _, _] = np.linalg.svd(multilinearalgebra.unfold(u, i))
                u = u[:, :ranks[i]]
                u = u.conj().T
                U.append(u)

            # Creating the core tensor
            S = multilinearalgebra.ten_mat_multiprod(tensor, np.array(U), *np.arange(0, tensor.ndim))

            # Convergence
            if multilinearalgebra.normalized_mean_square_error(S, S_init) > epsilon:
                S_init = S
                U_init = U

            else:
                for i in range(0, tensor.ndim):
                    U[i] = U[i].conj().T
                break

    k = np.array(U)[-2:]
    U = np.append(k, (np.array(U)[0:-2:]))

    return S, U, S.shape


def MLSKRF_3D(tensor, I, R):
    matrices = []
    for i in range(0, len(I)):
        a = np.zeros([I[i], R], dtype='complex')
        matrices.append(a)

    for i in range(0, R):
        tensor_aux = (tensor[:, i]).reshape(I)
        tensor_aux = np.moveaxis(tensor_aux, 1, 2)

        [S, U, _] = multilinearalgebra.HOSVD(tensor_aux)
        for j in range(0, len(I)):
            matrices[j][:, i] = (sqrt(S[0, 0, 0]) ** (2 / len(I))) * U[len(I) - j - 1][:, 0]

        tensor_approx = tl.tenalg.khatri_rao(matrices)
        nmse = normalized_mean_square_error(tensor, tensor_approx)

    return matrices, nmse


def MLSKRONF_3D(tensor, shapes, sizes):
    matrix_bar = extract_block_mlskronf(tensor, shapes)
    matrix_bar = vec(matrix_bar)
    tensor_bar = np.moveaxis(matrix_bar.reshape(sizes), len(shapes) - 2, len(shapes) - 1)
    [S, U, _] = multilinearalgebra.HOSVD(tensor_bar)

    matrices = []
    for i in range(0, tensor_bar.ndim):
        s = sqrt(S.flat[0]) ** (2 / len(shapes))
        u = s * U[::-1][i][:, 0].reshape(shapes[i], order='F')
        matrices.append(u)

    tensor_approx = tl.tenalg.kronecker(matrices)
    nmse = multilinearalgebra.normalized_mean_square_error(tensor, tensor_approx)

    return matrices, nmse


def CPD_ALS(tensor, rank, inter_max, delta):
    def unfold(tensor, n):

        if n == 0:
            tensor = np.moveaxis(tensor, n, 0)
            tensor_unfolding = tensor.reshape(tensor.shape[0], -1, order='F')

        else:
            tensor = np.moveaxis(tensor, n, 0)
            tensor_unfolding = tensor.reshape(tensor.shape[0], -1)

        return tensor_unfolding

    shapes = tensor.shape

    B = np.random.rand(shapes[2], rank)
    C = np.random.rand(shapes[0], rank)

    tensor_1mode = unfold(tensor, 1)
    tensor_2mode = unfold(tensor, 2)
    tensor_3mode = unfold(tensor, 0)

    error = np.zeros([inter_max])
    error[0] = 0
    for inter in range(1, inter_max):

        A = tensor_1mode @ np.linalg.pinv(((tl.tenalg.khatri_rao([C, B])).T))
        B = tensor_2mode @ np.linalg.pinv(((tl.tenalg.khatri_rao([C, A])).T))
        C = tensor_3mode @ np.linalg.pinv(((tl.tenalg.khatri_rao([B, A])).T))

        error[inter] = np.linalg.norm(tensor_1mode - A @ (tl.tenalg.khatri_rao([C, B]).T), 'fro') ** 2 / (
                np.linalg.norm(tensor_1mode, 'fro') ** 2)

        if abs(error[inter] - error[inter - 1]) <= delta:
            break

        else:
            continue

    return A, B, C, error[1:inter], inter - 1
