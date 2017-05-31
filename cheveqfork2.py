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

def Gauss(a, b0, b1):
    i_cur_stage = 0
    j_cur_stage = 0
    while j_cur_stage < 21 * 21:
        non_zero_ind = i_cur_stage
        while non_zero_ind < 21 * 21 - 1 and a[non_zero_ind, j_cur_stage] == 0:
            non_zero_ind += 1
        if a[non_zero_ind, j_cur_stage] == 0:
            j_cur_stage += 1
        else:
            for j in range(21 * 21):
                a[i_cur_stage, j], a[non_zero_ind, j] = a[non_zero_ind, j], a[i_cur_stage, j]
                b0[i_cur_stage, j], b0[non_zero_ind, j] = b0[non_zero_ind, j], b0[i_cur_stage, j]
                b1[i_cur_stage, j], b1[non_zero_ind, j] = b1[non_zero_ind, j], b1[i_cur_stage, j]
            for i in range(non_zero_ind + 1, 21 * 21):
                if a[i, j_cur_stage] != 0:
                    k = a[i, j_cur_stage] / a[i_cur_stage, j_cur_stage]
                    if k - int(k) != 0:
                        print("non int")
                    else:
                        k = int(k)
                    for j in range(21 * 21):
                        a[i, j] -= k * a[i_cur_stage, j]
            i_cur_stage += 1
            j_cur_stage += 1
    print(a[i_cur_stage - 1])
    return (a, b0, b1)



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
    #print(i)
    #print(w[i])
    q.append(w[i] * x[i])
    if not np.all(q[i] * q[i] * q[i] == np.identity(21, dtype=int)):
        print("=(",i)
    
for i1 in range(18):
    for j1 in range(18):
        for i2 in range(18):
            for j2 in range(18):
                for k in range(18):
                    if np.all(x[i1] * x[j1] * np.linalg.inv(x[i1]) * np.linalg.inv(x[j1]) == x[i2] * x[j2] * x[k]):
                        print(i1, "  , ", j1, " =  ", i2, " * ", j2, " * ", k)


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
counter = set()
for i in range(18):
    for j in range(18):
        sum_of_roots_i_j = roots[i] + roots[j]
        sum_of_roots_i_2j = roots[i] + 2 * roots[j]
        for ind, r in enumerate(roots):
            if np.array_equal(r, sum_of_roots_i_j):
                if np.all(x[i] * x[j] * np.linalg.inv(x[i]) * np.linalg.inv(x[j]) == x[ind]):
                    counter.add(i)
                    counter.add(j)
                    counter.add(ind)
                    #print(3, roots[i], roots[j], ind)
            for ind2, r2 in enumerate(roots):
                if np.array_equal(r2, sum_of_roots_i_2j):
                    if np.all(x[i] * x[j] * np.linalg.inv(x[i]) * np.linalg.inv(x[j]) == x[ind] * x[ind2]):
                        counter.add(i)
                        counter.add(j)
                        counter.add(ind)
                        counter.add(ind2)
                        #print(4, i, j, ind, ind2)
# print(counter)
if np.all(x[0] * x[2] * np.linalg.inv(x[0]) * np.linalg.inv(x[2]) == x[6]):
    print("win")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#linear system

var = 6
a = 0
b = 1
mc = []
base_factor = [0, 1, 6, 1, 6, 0, 6, 1, 6, 1]
if_inv = [False, False, False, True, True, True, False, False, True, True]

for i in range(10):
    left_factor = np.identity(21, dtype=int)
    right_factor = np.identity(21, dtype=int)
    #print(i)
    s = "L: "
    for l in range(i):
        if if_inv[l]:
            left_factor = left_factor * int_inv(x[base_factor[l]], 21, 21)
            s += str(base_factor[l]) + " -1, "
        else:
            left_factor = left_factor * x[base_factor[l]]
            s += str(base_factor[l]) + ", "
    #print(s)
    s = "R: "
    for r in range(i + 1, 10):
        if if_inv[r]:
            right_factor = right_factor * int_inv(x[base_factor[r]], 21, 21)
            s += str(base_factor[r]) + " -1, "
        else:
            right_factor = right_factor * x[base_factor[r]]
            s += str(base_factor[r]) + ", "
    #print(s)
    if if_inv[i]:
        left_factor = left_factor * int_inv(x[base_factor[i]], 21, 21)
        right_factor = int_inv(x[base_factor[i]], 21, 21) * right_factor
    
    mc.append(left_factor)
    mc.append(right_factor)

aa = np.zeros((21*21, 21*21), dtype=int)
bb0 = np.zeros((21*21, 21*21), dtype=int)
bb1 = np.zeros((21*21, 21*21), dtype=int)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_6
                aa[ii*21+jj, i*21+j] = mc[0][ii,i]*mc[0][j,jj] + mc[4][ii,i]*mc[5][j,jj] + mc[8][ii,i]*mc[9][j,jj] + mc[12][ii,i]*mc[13][j,jj] + mc[16][ii,i]*mc[17][j,jj]
                aa[ii*21+jj, i*21+j] = aa[ii * 21 + jj, i * 21 + j] % 2
                bb0[ii*21+jj, i*21+j] = (mc[0][ii,i]*mc[1][j,jj] + mc[10][ii,i]*mc[11][j,jj])
                bb1[ii*21+jj, i*21+j] = (mc[2][ii,i]*mc[3][j,jj] + mc[6][ii,i]*mc[7][j,jj] + mc[14][ii,i]*mc[15][j,jj] + mc[18][ii,i]*mc[19][j,jj])

# print(np.linalg.det(aa))


aa, bb0, bb1 = Gauss(aa, bb0, bb1)

for i in range(21*21):
    ss = str(i) + ": "
    flag = False
    for j in range(21*21):
        if aa[i, j] != 0:
            flag = True
            ss += str(j) + ", "
    if flag:
        #print(ss)
        pass


