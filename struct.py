roots = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
            (1, 1, 0), (-1, -1, 0), (0, 1, 1), (0, -1, -1), (0, 1, 2), (0, -1, -2),
            (1, 1, 1), (-1, -1, -1), (1, 1, 2), (-1, -1, -2), (1, 2, 2), (-1, -2, -2)]

def sum_root(a, b, c=(0, 0, 0)):
    return (a[0] + b[0] + c[0], a[1] + b[1] + c[1], a[2] + b[2] + c[2])

def minus_root(a):
    return (-a[0], -a[1], -a[2])

def sqr_root(a):    #go to basis
    eps = (a[0], a[1] - a[0], a[2] - a[1])
    return eps[0] * eps[0] + eps[1] * eps[1] + eps[2] * eps[2]

struct_const = [[ 0,  0,  1,  0,  0,  0,  0, -1, -1,  0, -1,  0,  0,  1,  0,  1,  0,  0],   #(1, 0, 0)
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
                [ 0,  0, -1,  0,  0,  0, -1,  0, -1,  0,  1,  0,  1,  0,  1,  0,  0,  0]]

for i in range(18):
    minus_i = roots.index(minus_root(roots[i]))
    for j in range(18):
        minus_j = roots.index(minus_root(roots[j]))
        if struct_const[i][j] != -struct_const[j][i]:
            print('bug 1 ', i, ' ', j)
        if struct_const[i][j] != -struct_const[minus_i][minus_j]:
            print('bug 2 ', i, ' ', j)

for i in range(18):
    for j in range(i + 1, 18):
        if sum_root(roots[i], roots[j]) in roots:
            a = roots[i]
            b = roots[j]
            c = minus_root(sum_root(a, b))
            k = roots.index(c)
            if sqr_root(a) * struct_const[i][j] != sqr_root(c) * struct_const[j][k]:
                print('bug 3 ', i, ' ', j)
            if sqr_root(b) * struct_const[j][k] != sqr_root(a) * struct_const[k][i]:
                print('bug 4 ', i, ' ', j)
        else:
            if struct_const[i][j]:
                print('bug 5 ', i, ' ', j)

for i in range(18):
    for j in range(i + 1, 18):
        for k in range(j + 1, 18):
            if sum_root(roots[i], roots[j], roots[k]) in roots:
                a = roots[i]
                b = roots[j]
                c = roots[k]
                d = minus_root(sum_root(a, b, c))
                t = roots.index(d)
                if sqr_root(sum_root(a, b)) and sqr_root(sum_root(b, c)) and sqr_root(sum_root(a, c)):
                    T1 = struct_const[i][j] * struct_const[k][t] / sqr_root(sum_root(a, b))
                    T2 = struct_const[j][k] * struct_const[i][t] / sqr_root(sum_root(b, c))
                    T3 = struct_const[k][i] * struct_const[j][t] / sqr_root(sum_root(c, a))
                    if T1 + T2 + T3:
                        print('bug 6 ', i, ' ', j, ' ', k )


'''i = 0
j = 2
if sum_root(roots[i], roots[j]) in roots:
    a = roots[i]
    b = roots[j]
    c = minus_root(sum_root(a, b)) # 1 * -1   != 2 * =((
    if sqr_root(a) * struct_const[roots.index(a)][roots.index(b)] != sqr_root(c) * struct_const[roots.index(b)][roots.index(c)]:
        print(c)'''
