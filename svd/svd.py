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
    for i in range (n+1) :
        e[i] = g
        s = 0
        l = i+1
        for j in range (i,m+1) :
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
            for j in range (l,n+1) :
                s = 0
                for k in range (i,m+1) :
                    s = s+u[k,i]*u[k,j]
                f = s/h
                for k in range (i,m+1) :
                    u[k,j] = u[k,j] + f*u[k,i]
        q[i] = g
        s = 0
        for j in range (l,n+1) :
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
            for j in range (l,n+1) :
                e[j] = u[i,j]/h
            for j in range (l,m+1) :
                s = 0
                for k in range (l,n+1) :
                    s = s+u[j,k]*u[i,k]
                for k in range (i,n+1) :
                    u[j,k] = u[j,k] + s*e[k]
        y = abs(q[i])+abs(e[i])
        if (y>x) :
            x = y
