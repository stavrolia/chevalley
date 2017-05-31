import numpy as np
from numpy import matrix

def scal_prod(a, b):
    a = a.tolist()[0]
    b = b.tolist()[0]
    eps_a = [a[0], a[1] - a[0], a[2] - a[1]]  #eps1-eps2, eps2-eps3, eps3
    eps_b = [b[0], b[1] - b[0], b[2] - b[1]]
    return eps_a[0] * eps_b[0] + eps_a[1] * eps_b[1] + eps_a[2] * eps_b[2]

def int_inv(a, m, n):
    inv_matr = np.linalg.inv(a)
    res = np.zeros((m, n), dtype=int)
    for i in range(m):
        for j in range(n):
            res[i, j] = int(inv_matr[i, j])
    return res

roots = matrix([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
            [1, 1, 0], [-1, -1, 0], [0, 1, 1], [0, -1, -1], [0, 1, 2], [0, -1, -2],
            [1, 1, 1], [-1, -1, -1], [1, 1, 2], [-1, -1, -2], [1, 2, 2], [-1, -2, -2]])

struct_const = matrix([[ 0,  0,  1,  0,  0,  0,  0, -1, -1,  0, -1,  0,  0,  1,  0,  1,  0,  0],   #(1, 0, 0)
                [ 0,  0,  0, -1,  0,  0,  1,  0,  0,  1,  0,  1, -1,  0, -1,  0,  0,  0],
                [-1,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0, -1,  0,  0,  1],   #(0, 1, 0)
                [ 0,  1,  0,  0,  0, -1, -1,  0,  1,  0,  0,  0,  0,  0,  0,  1, -1,  0],
                [ 0,  0, -1,  0,  0,  0,  1,  0,  2,  2,  0, -1,  2, -2,  0, -1,  0,  0],   #(0, 0, 1)
                [ 0,  0,  0,  1,  0,  0,  0, -1, -2, -2,  1,  0,  2, -2,  1,  0,  0,  0],
                [ 0, -1,  0,  1, -1,  0,  0,  0,  0,  0, -1,  0,  0,  1,  0,  0,  0,  1],   #(1, 1, 0)
                [ 1,  0, -1,  0,  0,  1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0, -1,  0],
                [ 1,  0,  0, -1, -2,  2,  0,  0,  0,  0,  0,  1, -2, -2,  0,  0,  0,  1],   #(0, 1, 1)
                [ 0, -1,  1,  0, -2,  2,  0,  0,  0,  0, -1,  0,  2,  2,  0,  0, -1,  0],
                [ 1,  0,  0,  0,  0, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0, -1],   #(0, 1, 2)
                [ 0, -1,  0,  0,  1,  0,  0, -1, -1,  0,  0,  0,  0,  0,  1,  0,  1,  0],
                [ 0,  1,  0,  0, -2, -2,  0,  1,  2, -2,  0,  0,  0,  0,  0,  1,  0, -1],   #(1, 1, 1)
                [-1,  0,  0,  0,  2,  2, -1,  0,  2, -2,  0,  0,  0,  0, -1,  0,  1,  0],
                [ 0,  1,  1,  0,  0, -1,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0, -1],   #(1, 1, 2)
                [-1,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1,  0],
                [ 0,  0,  0,  1,  0,  0,  0,  1,  0,  1,  0, -1,  0, -1,  0, -1,  0,  0],   #(1, 2, 2)
                [ 0,  0, -1,  0,  0,  0, -1,  0, -1,  0,  1,  0,  1,  0,  1,  0,  0,  0]])

adx = []
x = []
h = []
w = []
q = []

for i in range(3):
    h.append(np.zeros((21, 21), dtype=int))

for i in range(18):

    cur_adx = matrix(np.zeros((21,21), dtype=int))
    for j in range(18):
        sum_of_roots = roots[i] + roots[j]
        for ind, r in enumerate(roots):
            if np.array_equal(r, sum_of_roots):         #if sum_of_roots in roots:
                # print(roots[i], ' ', roots[j], ' ', sum_of_roots, ' ', ind)
                cur_adx[ind, j] = struct_const[i, j]

    minus_alpha = i + 1 - 2 * (i % 2)   # if i%2==0 then i + 1 else i - 1
    for j in range(3):
        cur_adx[i, 18 + j] = - 2 * scal_prod(roots[i], roots[2 * j]) // scal_prod(roots[2 * j], roots[2 * j])
        h[j][i, i] = 2 * scal_prod(roots[i], roots[2 * j]) // scal_prod(roots[2 * j], roots[2 * j])
        cur_adx[18 + j, minus_alpha] = (scal_prod(roots[2 * j], roots[2 * j]) * roots[i, j]) // scal_prod(roots[i], roots[i])

    if not np.all(cur_adx * cur_adx * cur_adx == np.zeros((21, 21), dtype=int)):
        print(i)
    adx.append(cur_adx)
    x.append(np.identity(21, dtype=int) + cur_adx + (cur_adx * cur_adx) // 2 )


for i in range(18):
    if i % 2 == 0:
        minus_alpha = i + 1
    else:
        minus_alpha = i - 1
    x_minus_alpha = np.identity(21, dtype=int) - adx[minus_alpha] + (adx[minus_alpha] * adx[minus_alpha]) // 2
    # print((adx[minus_alpha] * adx[minus_alpha]) / 2)
    w.append(x[i] * x_minus_alpha * x[i])
    print("w_", i, ":")
    print(w[i])
    q.append(w[i] * x[i])
    # if not np.all(q[i] * q[i] * q[i] == np.identity(21, dtype=int)):
    #     print("=(",i)
