import math

def SVD(m,n,withu,withv,eps,tol,a) :
    '''
    Kamus
    int : i,j,k,l,l1
    float : c,f,g,h,s,x,y,z
    array e
    '''
    q = [0 in range (n+1)]
    u = [[0 for j in range (n+1)] for i in range (m+1)]
    v = [[0 for j in range (n+1)] for i in range (n+1)]
    e = [0 for i in range (n+1)]
    for i in range (m+1) :
        for j in range (n+1) :
            u[i][j] = a[i,j]
    #Householder's reduction to bidiagonal form
    g = x = 0
    for i in range (n) :
        e[i] = g
        s = 0
        l = i+1
        for j in range (i,m) :
            s = s+(u[j,i])**2
        if s < tol :
            g = 0
        else :
            f = u[i,i]
            if (f<0) :
                g = math.sqrt(s)
            else :
                g = -math.sqrt(s)
            h = f*g-s
            u[i,i] = f-g
            for j in range (l,n) :
                s = 0
                for k in range (i,m) :
                    s = s+u[k,i]*u[k,j]
                f = s/h
                for k in range (i,m) :
                    u[k,j] = u[k,j] + f*u[k,i]
        q[i] = g
        s = 0
        for j in range (l,n) :
            s = s+(u[i,j])**2
        if (s<tol) :
            g = 0
        else :
            f = u[i,i+1]
            if (f<0) :
                g = math.sqrt(s)
            else :
                g = -math.sqrt(s)
            h = f*g-s
            u[i,i+1] = f-g
            for j in range (l,n) :
                e[j] = u[i,j]/h
            for j in range (l,m) :
                s = 0
                for k in range (l,n) :
                    s = s+u[j,k]*u[i,k]
                for k in range (i,n) :
                    u[j,k] = u[j,k] + s*e[k]
        y = abs(q[i])+abs(e[i])
        if (y>x) :
            x = y

# COMMENT 3
# accumulation of right-hand transformations

    for i in range(n-1, -1, -1):
        if (g != 0):
            h = u[i][i+1] * g
            for j in range(l, n):
                v[j][i] = u[i][j]/h
            for j in range(l, n):
                s = 0
                for k in range(l, n):
                    s = s + u[i][k] * v[k][j]
                for k in range(l, n):
                    v[k][j] = v[k][j] + s * v[k][i]
        for j in range(l, n):
            v[j][i] = 0
            v[i][j] = v[j][i] 
        v[i][i] = 1
        g = e[i]
        l = i    

    # COMMENT 4
    # accumulation of left-hand transformations; 

    for i in range(n-1, -1, -1):
        l = i + 1
        g = q[i]
        for j in range(l, n):
            u[i][j] = 0
        if (g != 0):
            h = u[i][i] * g 
            for j in range(l, n):
                s = 0
                for k in range(l, m):
                    s = s + u[k][i] * u[k][j]
                f = s/h
                for k in range(i, m):
                    u[k][j] = u[k][j] + f * u[k][i]
            for j in range(i, m):
                u[j][i] = u[j][i]/g
        else :
            for j in range(i, m):
                u[j][i] = 0
            u[i][i] = u[i][i] + 1



