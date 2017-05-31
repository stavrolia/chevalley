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
            if inv_matr[i, j] - res[i, j] != 0:
                print("=(")
    return res

def count_const(a, b):
    res = []
    res.append(np.identity(21, dtype=int))
    res.append(x[b] * int_inv(x[a], 21, 21) * int_inv(x[b], 21, 21))
    res.append(x[a])
    res.append(int_inv(x[a], 21, 21) * int_inv(x[b], 21, 21))
    res.append(x[a] * x[b] * int_inv(x[a], 21, 21))
    res.append(int_inv(x[a], 21, 21) * int_inv(x[b], 21, 21))
    res.append(x[a] * x[b] * int_inv(x[a], 21, 21) * int_inv(x[b], 21, 21))
    res.append(int_inv(x[b], 21, 21))
    return res

def sh(a):
    shifts = [False, False, 0*21*21, 21*21, False, False, 2*21*21, 3*21*21, 8*21*21, 9*21*21, 4*21*21, 5*21*21, 10*21*21, 11*21*21, 6*21*21, 7*21*21]
    return shifts[a]

def Gauss(a, b0, b1, b4, b5, b16, b17): #aa * x = bb0 + bb1 + bb4 + bb5 + bb16 + bb17
    i_cur_stage = 0
    j_cur_stage = 0
    n, m = np.shape(a)

    while j_cur_stage < m:
        non_zero_ind = i_cur_stage
        while (non_zero_ind < n - 1) and (a[non_zero_ind, j_cur_stage] == 0):
            non_zero_ind += 1
        if a[non_zero_ind, j_cur_stage] == 0:
            j_cur_stage += 1
        else:
            for j in range(m):
                a[i_cur_stage, j], a[non_zero_ind, j] = a[non_zero_ind, j], a[i_cur_stage, j]
            for j in range(21 * 21):
                b0[i_cur_stage, j], b0[non_zero_ind, j] = b0[non_zero_ind, j], b0[i_cur_stage, j]
                b1[i_cur_stage, j], b1[non_zero_ind, j] = b1[non_zero_ind, j], b1[i_cur_stage, j]
                b4[i_cur_stage, j], b4[non_zero_ind, j] = b4[non_zero_ind, j], b4[i_cur_stage, j]
                b5[i_cur_stage, j], b5[non_zero_ind, j] = b5[non_zero_ind, j], b5[i_cur_stage, j]
                b16[i_cur_stage, j], b16[non_zero_ind, j] = b16[non_zero_ind, j], b16[i_cur_stage, j]
                b17[i_cur_stage, j], b17[non_zero_ind, j] = b17[non_zero_ind, j], b17[i_cur_stage, j]
            for i in range(non_zero_ind + 1, n):
                if a[i, j_cur_stage] != 0:
                    for j in range(m):
                        a[i, j] -= a[i_cur_stage, j]
                        a[i, j] = a[i, j] % 2
            i_cur_stage += 1
            j_cur_stage += 1
    return (a, b0, b1, b4, b5, b16, b17)

def second_Gauss(a, n): 
    i_cur_stage = 0
    j_cur_stage = 0
    m = 6 * 21 * 21

    while j_cur_stage < m:
        non_zero_ind = i_cur_stage
        while (non_zero_ind < n - 1) and (a[non_zero_ind, j_cur_stage] == 0):
            non_zero_ind += 1
        if a[non_zero_ind, j_cur_stage] == 0:
            j_cur_stage += 1
        else:
            for j in range(m):
                a[i_cur_stage, j], a[non_zero_ind, j] = a[non_zero_ind, j], a[i_cur_stage, j]
            for i in range(non_zero_ind + 1, n):
                if a[i, j_cur_stage] != 0:
                    for j in range(m):
                        a[i, j] -= a[i_cur_stage, j]
                        a[i, j] = a[i, j] % 2
            i_cur_stage += 1
            j_cur_stage += 1
    return a


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

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#linear system



aa = np.zeros((60*21*21, 12*21*21), dtype=int)  #aa * x = bb0 + bb1 + bb4 + bb5 + bb16 + bb17
bb0 = np.zeros((60*21*21, 21*21), dtype=int)
bb1 = np.zeros((60*21*21, 21*21), dtype=int)
bb4 = np.zeros((60*21*21, 21*21), dtype=int)
bb5 = np.zeros((60*21*21, 21*21), dtype=int)
bb16 = np.zeros((60*21*21, 21*21), dtype=int)
bb17 = np.zeros((60*21*21, 21*21), dtype=int)

#first equation, shift i == 0
first = count_const(0, 2)    #const for equation from index of X in commutator
second = count_const(0, 15)
third = count_const(16, 3)
fourth = count_const(16, 7)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[ii * 21 + jj, sh(6) + i * 21 + j] = first[0][ii,i] * first[0][j,jj] % 2 #X_6
                aa[6*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = second[0][ii,i] * second[0][j,jj] % 2 #X_11
                aa[12*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = third[0][ii,i] * third[0][j,jj] % 2 #X_14
                aa[18*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = fourth[0][ii,i] * fourth[0][j,jj] % 2 #X_10
                #second argument of commutator
                aa[ii * 21 + jj, sh(2) + i * 21 + j] = (- first[2][ii, i] * first[3][j,jj] - first[6][ii,i] * first[7][j,jj]) % 2 #X_2
                aa[6*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = (- second[2][ii,i] * second[3][j,jj] - second[6][ii,i] * second[7][j,jj]) % 2 #X_15
                aa[12*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = (- third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2 #X_3
                aa[18*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = (- fourth[2][ii,i] * fourth[3][j,jj] - fourth[6][ii,i] * fourth[7][j,jj]) % 2 #X_7
                #first argument
                bb0[ii * 21 + jj, i * 21 + j] = (first[0][ii,i] * first[1][j,jj] + first[4][ii,i] * first[5][j,jj]) % 2
                bb0[6*21*21 + ii * 21 + jj, i * 21 + j] = (second[0][ii,i] * second[1][j,jj] + second[4][ii,i] * second[5][j,jj]) % 2
                bb16[12*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                bb16[18*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[0][ii,i] * fourth[1][j,jj] + fourth[4][ii,i] * fourth[5][j,jj]) % 2

#second equation, shift i == 21*21
first = count_const(1, 6)
second = count_const(1, 11)
third = count_const(17, 14)
fourth = count_const(17, 10)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = first[0][ii,i] * first[0][j,jj] % 2
                aa[7*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = second[0][ii,i] * second[0][j,jj] % 2
                aa[13*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = third[0][ii,i] * third[0][j,jj] % 2
                aa[19*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = fourth[0][ii,i] * fourth[0][j,jj] % 2
                #second argument
                aa[21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (- first[2][ii, i] * first[3][j,jj] - first[6][ii,i] * first[7][j,jj]) % 2
                aa[7*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (- second[2][ii, i] * second[3][j,jj] - second[6][ii,i] * second[7][j,jj]) % 2
                aa[13*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (- third[2][ii, i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[19*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = (- fourth[2][ii, i] * fourth[3][j,jj] - fourth[6][ii,i] * fourth[7][j,jj]) % 2
                #first argument
                bb1[21*21 + ii * 21 + jj, i * 21 + j] = (first[0][ii,i] * first[1][j,jj] + first[4][ii,i] * first[5][j,jj]) % 2
                bb1[7*21*21 + ii * 21 + jj, i * 21 + j] = (second[0][ii,i] * second[1][j,jj] + second[4][ii,i] * second[5][j,jj]) % 2
                bb17[13*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                bb17[19*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[0][ii,i] * fourth[1][j,jj] + fourth[4][ii,i] * fourth[5][j,jj]) % 2

#third equation, shift i = 2*21*21
first = count_const(3, 1)
second = count_const(14, 1)
third = count_const(2, 17)
fourth = count_const(6, 17)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[2*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = first[0][ii,i] * first[0][j,jj] % 2
                aa[8*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = second[0][ii,i] * second[0][j,jj] % 2
                aa[14*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = third[0][ii,i] * third[0][j,jj] % 2
                aa[20*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = fourth[0][ii,i] * fourth[0][j,jj] % 2
                #res
                aa[2*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = (- first[0][ii, i] * first[1][j,jj] - first[4][ii,i] * first[5][j,jj]) % 2
                aa[8*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (- second[0][ii, i] * second[1][j,jj] - second[4][ii,i] * second[5][j,jj]) % 2
                aa[14*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = (- third[0][ii, i] * third[1][j,jj] - third[4][ii,i] * third[5][j,jj]) % 2
                aa[20*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (- fourth[0][ii, i] * fourth[1][j,jj] - fourth[4][ii,i] * fourth[5][j,jj]) % 2
                #second arg
                bb1[2*21*21 + ii * 21 + jj, i * 21 + j] = (first[2][ii,i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                bb1[8*21*21 + ii * 21 + jj, i * 21 + j] = (second[2][ii,i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2
                bb17[14*21*21 + ii * 21 + jj, i * 21 + j] = (third[2][ii,i] * third[3][j,jj] + third[6][ii,i] * third[7][j,jj]) % 2
                bb17[20*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[2][ii,i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2

#fourth eq, shift i = 3*21*21
first = count_const(7, 0)
second = count_const(10, 0)
third = count_const(15, 16)
fourth = count_const(11, 16)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[3*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = first[0][ii,i] * first[0][j,jj] % 2
                aa[9*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = second[0][ii,i] * second[0][j,jj] % 2
                aa[15*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = third[0][ii,i] * third[0][j,jj] % 2
                aa[21*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = fourth[0][ii,i] * fourth[0][j,jj] % 2
                #res
                aa[3*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = (- first[0][ii, i] * first[1][j,jj] - first[4][ii,i] * first[5][j,jj]) % 2
                aa[9*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = (- second[0][ii, i] * second[1][j,jj] - second[4][ii,i] * second[5][j,jj]) % 2
                aa[15*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = (- third[0][ii, i] * third[1][j,jj] - third[4][ii,i] * third[5][j,jj]) % 2
                aa[21*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (- fourth[0][ii, i] * fourth[1][j,jj] - fourth[4][ii,i] * fourth[5][j,jj]) % 2
                #second arg
                bb0[3*21*21 + ii * 21 + jj, i * 21 + j] = (first[2][ii,i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                bb0[9*21*21 + ii * 21 + jj, i * 21 + j] = (second[2][ii,i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2
                bb16[15*21*21 + ii * 21 + jj, i * 21 + j] = (third[2][ii,i] * third[3][j,jj] + third[6][ii,i] * third[7][j,jj]) % 2
                bb16[21*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[2][ii,i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2

#fifth eq, shift i = 4*21*21
first = count_const(2, 7)
second = count_const(15, 10)
third = count_const(3, 15)
fourth = count_const(7, 11)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[4*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = (first[0][ii,i] * first[1][j,jj] + first[4][ii,i] * first[5][j,jj]) % 2
                aa[10*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = (second[0][ii,i] * second[1][j,jj] + second[4][ii,i] * second[5][j,jj]) % 2
                aa[16*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                aa[22*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = (fourth[0][ii,i] * fourth[1][j,jj] + fourth[4][ii,i] * fourth[5][j,jj]) % 2
                #second arg
                aa[4*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = (first[2][ii, i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                aa[10*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = (second[2][ii, i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2
                aa[16*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = (third[2][ii, i] * third[3][j,jj] + third[6][ii,i] * third[7][j,jj]) % 2
                aa[22*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (fourth[2][ii, i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2
                # res
                bb1[4*21*21 + ii * 21 + jj, i * 21 + j] = first[0][ii,i] * first[0][j,jj] % 2
                bb1[10*21*21 + ii * 21 + jj, i * 21 + j] = second[0][ii,i] * second[0][j,jj] % 2
                bb17[16*21*21 + ii * 21 + jj, i * 21 + j] = third[0][ii,i] * third[0][j,jj] % 2
                bb17[22*21*21 + ii * 21 + jj, i * 21 + j] = fourth[0][ii,i] * fourth[0][j,jj] % 2

#sixth eq, shift i = 5*21*21
first = count_const(6, 3)
second = count_const(11, 14)
third = count_const(14, 2)
fourth = count_const(10, 6)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[5*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (first[0][ii,i] * first[1][j,jj] + first[4][ii,i] * first[5][j,jj]) % 2
                aa[11*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (second[0][ii,i] * second[1][j,jj] + second[4][ii,i] * second[5][j,jj]) % 2
                aa[17*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                aa[23*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = (fourth[0][ii,i] * fourth[1][j,jj] + fourth[4][ii,i] * fourth[5][j,jj]) % 2
                #second arg
                aa[5*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = (first[2][ii, i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                aa[11*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (second[2][ii, i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2
                aa[17*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = (third[2][ii, i] * third[3][j,jj] + third[6][ii,i] * third[7][j,jj]) % 2
                aa[23*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (fourth[2][ii, i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2
                # res
                bb0[5*21*21 + ii * 21 + jj, i * 21 + j] = first[0][ii,i] * first[0][j,jj] % 2
                bb0[11*21*21 + ii * 21 + jj, i * 21 + j] = second[0][ii,i] * second[0][j,jj] % 2
                bb16[17*21*21 + ii * 21 + jj, i * 21 + j] = third[0][ii,i] * third[0][j,jj] % 2
                bb16[23*21*21 + ii * 21 + jj, i * 21 + j] = fourth[0][ii,i] * fourth[0][j,jj] % 2

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
first = count_const(8, 11)
second = count_const(9, 2)
third = count_const(6, 13)
fourth = count_const(14, 3)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[24*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = first[0][ii,i] * x[5][j,jj] % 2
                aa[25*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = x[5][ii,i] * first[0][j,jj] % 2
                aa[26*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = second[0][ii,i] * x[5][j,jj] % 2
                aa[27*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = x[5][ii,i] * second[0][j,jj] % 2
                aa[28*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = third[0][ii,i] * x[5][j,jj] % 2
                aa[29*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = x[5][ii,i] * third[0][j,jj] % 2
                aa[30*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = fourth[0][ii,i] * x[4][j,jj] % 2
                aa[31*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = x[4][ii,i] * fourth[0][j,jj] % 2

                aa[24*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = (-first[0][ii,i] * first[1][j,jj] - first[4][ii,i] * first[5][j,jj]) % 2
                aa[25*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = (-first[0][ii,i] * first[1][j,jj] - first[4][ii,i] * first[5][j,jj]) % 2
                aa[26*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (-second[0][ii,i] * second[1][j,jj] - second[4][ii,i] * second[5][j,jj]) % 2
                aa[27*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (-second[0][ii,i] * second[1][j,jj] - second[4][ii,i] * second[5][j,jj]) % 2
                aa[28*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (-third[0][ii,i] * third[1][j,jj] - third[4][ii,i] * third[5][j,jj]) % 2
                aa[29*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (-third[0][ii,i] * third[1][j,jj] - third[4][ii,i] * third[5][j,jj]) % 2
                aa[30*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (-fourth[0][ii,i] * fourth[1][j,jj] - fourth[4][ii,i] * fourth[5][j,jj]) % 2
                aa[31*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (-fourth[0][ii,i] * fourth[1][j,jj] - fourth[4][ii,i] * fourth[5][j,jj]) % 2

                aa[24*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (-first[2][ii,i] * first[3][j,jj] - first[6][ii,i] * first[7][j,jj]) % 2
                aa[25*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (-first[2][ii,i] * first[3][j,jj] - first[6][ii,i] * first[7][j,jj]) % 2
                aa[26*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = (-second[2][ii,i] * second[3][j,jj] - second[6][ii,i] * second[7][j,jj]) % 2
                aa[27*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = (-second[2][ii,i] * second[3][j,jj] - second[6][ii,i] * second[7][j,jj]) % 2
                aa[28*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (-third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[29*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (-third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[30*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (-fourth[2][ii,i] * fourth[3][j,jj] - fourth[6][ii,i] * fourth[7][j,jj]) % 2
                aa[31*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (-fourth[2][ii,i] * fourth[3][j,jj] - fourth[6][ii,i] * fourth[7][j,jj]) % 2

                bb5[24*21*21 + ii * 21 + jj, i * 21 + j] = (-x[2][ii, i] * first[0][j, jj]) % 2
                bb5[25*21*21 + ii * 21 + jj, i * 21 + j] = (-first[0][ii, i] * x[2][j, jj]) % 2
                bb5[26*21*21 + ii * 21 + jj, i * 21 + j] = (-x[11][ii, i] * second[0][j, jj]) % 2
                bb5[27*21*21 + ii * 21 + jj, i * 21 + j] = (-second[0][ii, i] * x[11][j, jj]) % 2
                bb5[28*21*21 + ii * 21 + jj, i * 21 + j] = (-x[15][ii, i] * third[0][j, jj]) % 2
                bb5[29*21*21 + ii * 21 + jj, i * 21 + j] = (-third[0][ii, i] * x[15][j, jj]) % 2
                bb4[30*21*21 + ii * 21 + jj, i * 21 + j] = (-x[7][ii, i] * fourth[0][j, jj]) % 2
                bb4[31*21*21 + ii * 21 + jj, i * 21 + j] = (-fourth[0][ii, i] * x[7][j, jj]) % 2

first = count_const(2, 4)
second = count_const(11, 4)
third = count_const(5, 14)
fourth = count_const(4, 6)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[32*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = first[0][ii,i] * x[10][j,jj] % 2
                aa[33*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = x[10][ii,i] * first[0][j,jj] % 2
                aa[34*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = second[0][ii,i] * x[10][j,jj] % 2
                aa[35*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = x[10][ii,i] * second[0][j,jj] % 2

                aa[32*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = x[8][ii, i] * first[0][j, jj] % 2
                aa[33*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = first[0][ii, i] * x[8][j, jj] % 2
                aa[34*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = x[8][ii, i] * second[0][j, jj] % 2
                aa[35*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = second[0][ii, i] * x[8][j, jj] % 2

                aa[32*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = (-first[0][ii,i] * first[1][j,jj] - first[4][ii,i] * first[5][j,jj]) % 2
                aa[33*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = (-first[0][ii,i] * first[1][j,jj] - first[4][ii,i] * first[5][j,jj]) % 2
                aa[34*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (-second[0][ii,i] * second[1][j,jj] - second[4][ii,i] * second[5][j,jj]) % 2
                aa[35*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (-second[0][ii,i] * second[1][j,jj] - second[4][ii,i] * second[5][j,jj]) % 2
                                
                bb4[32*21*21 + ii * 21 + jj, i * 21 + j] = (first[2][ii,i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                bb4[33*21*21 + ii * 21 + jj, i * 21 + j] = (first[2][ii,i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                bb4[34*21*21 + ii * 21 + jj, i * 21 + j] = (second[2][ii,i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2
                bb4[35*21*21 + ii * 21 + jj, i * 21 + j] = (second[2][ii,i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2


                aa[36*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = third[0][ii,i] * x[6][j,jj] % 2
                aa[37*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = x[6][ii,i] * third[0][j,jj] % 2
                aa[38*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = fourth[0][ii,i] * x[6][j,jj] % 2
                aa[39*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = x[6][ii,i] * fourth[0][j,jj] % 2

                aa[36*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = x[12][ii,i] * third[0][j,jj] % 2
                aa[37*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = third[0][ii,i] * x[12][j,jj] % 2
                aa[38*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = x[12][ii,i] * fourth[0][j,jj] % 2
                aa[39*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = fourth[0][ii,i] * x[12][j,jj] % 2

                aa[36*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (-third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[37*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (-third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[38*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (-fourth[2][ii,i] * fourth[3][j,jj] - fourth[6][ii,i] * fourth[7][j,jj]) % 2
                aa[39*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (-fourth[2][ii,i] * fourth[3][j,jj] - fourth[6][ii,i] * fourth[7][j,jj]) % 2

                bb5[36*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                bb5[37*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                bb4[38*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[0][ii,i] * fourth[1][j,jj] + fourth[4][ii,i] * fourth[5][j,jj]) % 2
                bb4[39*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[0][ii,i] * fourth[1][j,jj] + fourth[4][ii,i] * fourth[5][j,jj]) % 2

first = count_const(4, 8)
second = count_const(8, 5)
third = count_const(4, 9)
fourth = count_const(9, 5)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[40*21*21 + ii * 21 + jj, sh(10) + i * 21 + j] = (first[0][ii,i] * x[10][j,jj] + x[10][ii,i] * first[0][j,jj]) % 2
                aa[41*21*21 + ii * 21 + jj, sh(2) + i * 21 + j] = (second[0][ii,i] * x[2][j,jj] + x[2][ii,i] * second[0][j,jj]) % 2
                aa[42*21*21 + ii * 21 + jj, sh(3) + i * 21 + j] = (third[0][ii,i] * x[3][j,jj] + x[3][ii,i] * third[0][j,jj]) % 2
                aa[43*21*21 + ii * 21 + jj, sh(11) + i * 21 + j] = (fourth[0][ii,i] * x[11][j,jj] + x[11][ii,i] * fourth[0][j,jj]) % 2

                aa[40*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = (-first[2][ii,i] * first[3][j,jj] - first[6][ii,i] * first[7][j,jj]) % 2
                aa[41*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = (-second[0][ii,i] * second[1][j,jj] - second[4][ii,i] * second[5][j,jj]) % 2
                aa[42*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (-third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[43*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (-fourth[0][ii,i] * fourth[1][j,jj] - fourth[4][ii,i] * fourth[5][j,jj]) % 2

                bb4[40*21*21 + ii * 21 + jj, i * 21 + j] = (first[0][ii,i] * first[1][j,jj] + first[4][ii,i] * first[5][j,jj]) % 2
                bb5[41*21*21 + ii * 21 + jj, i * 21 + j] = (second[2][ii,i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2
                bb4[42*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                bb5[43*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[2][ii,i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2

first = count_const(5, 12)
second = count_const(13, 5)
third = count_const(4, 12)
fourth = count_const(13, 4)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[44*21*21 + ii * 21 + jj, sh(6) + i * 21 + j] = (first[0][ii,i] * x[6][j,jj] + x[6][ii,i] * first[0][j,jj]) % 2
                aa[45*21*21 + ii * 21 + jj, sh(15) + i * 21 + j] = (second[0][ii,i] * x[15][j,jj] + x[15][ii,i] * second[0][j,jj]) % 2
                aa[46*21*21 + ii * 21 + jj, sh(14) + i * 21 + j] = (third[0][ii,i] * x[14][j,jj] + x[14][ii,i] * third[0][j,jj]) % 2
                aa[47*21*21 + ii * 21 + jj, sh(7) + i * 21 + j] = (fourth[0][ii,i] * x[7][j,jj] + x[7][ii,i] * fourth[0][j,jj]) % 2

                aa[44*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = (-first[2][ii,i] * first[3][j,jj] - first[6][ii,i] * first[7][j,jj]) % 2
                aa[45*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (-second[0][ii,i] * second[1][j,jj] - second[4][ii,i] * second[5][j,jj]) % 2
                aa[46*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = (-third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[47*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (-fourth[0][ii,i] * fourth[1][j,jj] - fourth[4][ii,i] * fourth[5][j,jj]) % 2

                bb5[44*21*21 + ii * 21 + jj, i * 21 + j] = (first[0][ii,i] * first[1][j,jj] + first[4][ii,i] * first[5][j,jj]) % 2
                bb5[45*21*21 + ii * 21 + jj, i * 21 + j] = (second[2][ii,i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2
                bb4[46*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                bb4[47*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[2][ii,i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2

first = count_const(13, 16)
second = count_const(16, 9)
third = count_const(1, 9)
fourth = count_const(12, 1)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[48*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (-first[0][ii,i] * first[1][j,jj] - first[4][ii,i] * first[5][j,jj]) % 2
                aa[49*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (-first[0][ii,i] * first[1][j,jj] - first[4][ii,i] * first[5][j,jj]) % 2
                aa[50*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (-second[2][ii,i] * second[3][j,jj] - second[6][ii,i] * second[7][j,jj]) % 2
                aa[51*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (-second[2][ii,i] * second[3][j,jj] - second[6][ii,i] * second[7][j,jj]) % 2
                aa[52*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (-third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[53*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (-third[2][ii,i] * third[3][j,jj] - third[6][ii,i] * third[7][j,jj]) % 2
                aa[54*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = (-fourth[0][ii,i] * fourth[1][j,jj] - fourth[4][ii,i] * fourth[5][j,jj]) % 2
                aa[55*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = (-fourth[0][ii,i] * fourth[1][j,jj] - fourth[4][ii,i] * fourth[5][j,jj]) % 2

                aa[48*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = x[1][ii,i] * first[0][j,jj] % 2
                aa[49*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = first[0][ii,i] * x[1][j,jj] % 2
                aa[50*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = x[0][ii,i] * second[0][j,jj] % 2
                aa[51*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = second[0][ii,i] * x[0][j,jj] % 2
                aa[52*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = x[17][ii,i] * third[0][j,jj] % 2
                aa[53*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = third[0][ii,i] * x[17][j,jj] % 2
                aa[54*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = fourth[0][ii,i] * x[16][j,jj] % 2
                aa[55*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = x[16][ii,i] * fourth[0][j,jj] % 2

                bb1[48*21*21 + ii * 21 + jj, i * 21 + j] = -first[0][ii,i] * x[8][j,jj] % 2
                bb1[49*21*21 + ii * 21 + jj, i * 21 + j] = -x[8][ii,i] * first[0][j,jj] % 2
                bb0[50*21*21 + ii * 21 + jj, i * 21 + j] = -second[0][ii,i] * x[12][j,jj] % 2
                bb0[51*21*21 + ii * 21 + jj, i * 21 + j] = -x[12][ii,i] * second[0][j,jj] % 2
                bb17[52*21*21 + ii * 21 + jj, i * 21 + j] = -third[0][ii,i] * x[13][j,jj] % 2
                bb17[53*21*21 + ii * 21 + jj, i * 21 + j] = -x[13][ii,i] * third[0][j,jj] % 2
                bb16[54*21*21 + ii * 21 + jj, i * 21 + j] = -x[8][ii,i] * fourth[0][j,jj] % 2
                bb16[55*21*21 + ii * 21 + jj, i * 21 + j] = -fourth[0][ii,i] * x[8][j,jj] % 2

                bb16[48*21*21 + ii * 21 + jj, i * 21 + j] = (first[2][ii,i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                bb16[49*21*21 + ii * 21 + jj, i * 21 + j] = (first[2][ii,i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                bb16[50*21*21 + ii * 21 + jj, i * 21 + j] = (second[0][ii,i] * second[1][j,jj] + second[4][ii,i] * second[5][j,jj]) % 2
                bb16[51*21*21 + ii * 21 + jj, i * 21 + j] = (second[0][ii,i] * second[1][j,jj] + second[4][ii,i] * second[5][j,jj]) % 2
                bb1[52*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                bb1[53*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                bb1[54*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[2][ii,i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2
                bb1[55*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[2][ii,i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2

first = count_const(13, 8)
second = count_const(9, 12)
third = count_const(9, 13)
fourth = count_const(12, 8)

for ii in range(21):
    for jj in range(21):     # [ii, jj] - number of cell in big matrix
        for i in range(21):
            for j in range(21):     # [i, j] - index of X_k
                aa[56*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (first[0][ii,i] * first[1][j,jj] + first[4][ii,i] * first[5][j,jj]) % 2
                aa[57*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (second[0][ii,i] * second[1][j,jj] + second[4][ii,i] * second[5][j,jj]) % 2
                aa[58*21*21 + ii * 21 + jj, sh(9) + i * 21 + j] = (third[0][ii,i] * third[1][j,jj] + third[4][ii,i] * third[5][j,jj]) % 2
                aa[59*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = (fourth[0][ii,i] * fourth[1][j,jj] + fourth[4][ii,i] * fourth[5][j,jj]) % 2

                aa[56*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = (first[2][ii,i] * first[3][j,jj] + first[6][ii,i] * first[7][j,jj]) % 2
                aa[57*21*21 + ii * 21 + jj, sh(12) + i * 21 + j] = (second[2][ii,i] * second[3][j,jj] + second[6][ii,i] * second[7][j,jj]) % 2
                aa[58*21*21 + ii * 21 + jj, sh(13) + i * 21 + j] = (third[2][ii,i] * third[3][j,jj] + third[6][ii,i] * third[7][j,jj]) % 2
                aa[59*21*21 + ii * 21 + jj, sh(8) + i * 21 + j] = (fourth[2][ii,i] * fourth[3][j,jj] + fourth[6][ii,i] * fourth[7][j,jj]) % 2

                bb1[56*21*21 + ii * 21 + jj, i * 21 + j] = (first[0][ii,i] * x[1][j,jj] + x[1][ii,i] * first[0][j,jj]) % 2
                bb0[57*21*21 + ii * 21 + jj, i * 21 + j] = (second[0][ii,i] * x[0][j,jj] + x[0][ii,i] * second[0][j,jj]) % 2
                bb17[58*21*21 + ii * 21 + jj, i * 21 + j] = (third[0][ii,i] * x[17][j,jj] + x[17][ii,i] * third[0][j,jj]) % 2
                bb16[59*21*21 + ii * 21 + jj, i * 21 + j] = (fourth[0][ii,i] * x[16][j,jj] + x[16][ii,i] * fourth[0][j,jj]) % 2

aa, bb0, bb1, bb4, bb5, bb16, bb17 = Gauss(aa, bb0, bb1, bb4, bb5, bb16, bb17)

A = np.zeros((60*21*21, 6*21*21))
n_A = 0

for i in range(np.shape(aa)[0]):
    flag_a = False
    flag_b = False

    for j in range(np.shape(aa)[1]):
        if aa[i, j] != 0:
            flag_a = True
    for j in range(21 * 21):
        if bb0[i, j] or bb1[i, j] or bb4[i, j] or bb5[i,j] or bb16[i, j] or bb17[i, j]:
            flag_b = True
    if flag_b and not flag_a:
        for j in range(21 * 21):
            A[n_A, j] = bb0[i, j]
            A[n_A, 21*21 + j] = bb1[i, j]
            A[n_A, 2*21*21 + j] = bb4[i, j]
            A[n_A, 3*21*21 + j] = bb5[i, j]
            A[n_A, 4*21*21 + j] = bb16[i, j]
            A[n_A, 5*21*21 + j] = bb17[i, j]
        n_A += 1

A = second_Gauss(A, n_A)

f = open("chev2.txt", "w")
counter = 0

for i in range(n_A):
    flag = False
    s = str(i) + ": "
    for j in range(6 * 21 * 21):
        if A[i, j]:
            flag = True
            s = s + str(j) + ", "

    if flag:
        s = s + "\n"
        counter += 1
        f.write(s)
    
f.close()
print(counter)


