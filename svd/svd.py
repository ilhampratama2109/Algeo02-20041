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

    # COMMENT 5
    # diagonalization of tile bidiagonal form;

    eps = eps*x
    for k in range(n-1, -1, -1):
        for l in range(k-1, -1, -1):  #tes f splitting
            if (abs(e[l]) <= eps): #masuk ke test f konvergen
                z = float(q[k]) 
                if (l == k):    #masuk ke konvergen
                    if (z < 0):
                        q[k] = -z
                        for j in range(1, n):
                            v[j,k] = -v[j,k]
                
                x = q[l]
                y = q[k-1]
                g = e[k-1]
                h = e[k]
                f = ((y-x)*(y+z)*(g-h)*(g+h))/(2*h*y)
                g = math.sqrt(f*f+1)
                if(f < 0):
                    f = ((x-z)*(x+z)+h*(y/(f-g)-h))/x
                else:
                    f = ((x-z)*(x+z)+h*(y/(f+g)-h))/x
                
                c = s = 1
                for i in range(l+1,k):
                    g = e[i]
                    y = q[i]
                    h = s*g
                    g = c*g
                    e[i-1] = z = math.sqrt(f*f+h*h)
                    c = f/z
                    s = h/z
                    f = x*c+g*s
                    g = -x*s+g*c
                    h = y*s
                    y = y*c
                    for j in range(n+1):
                        x = v[j,i-1]
                        z = v[j,i]
                        v[j,i-1] = x*c+z*s
                        v[j,i] = -x*s+z*c
                    q[i-1] = z = math.sqrt(f*f+h*h)
                    c = f/z
                    s = h/z
                    f = c*g+s*y
                    x = -s*g+c*y
                    for j in range(m+1):
                        y = u[j, i-1]
                        z = u[j,i]
                        u[j,i-1] = y*c+z*s
                        u[j,i] = -y*s+z*c
                e[l] = 0
                e[k] = f
                q[k] = x
            
            if(abs(q[l-1]) <= eps): #masuk ke cancellation
                c = 0 
                s = 1
                l1 = l-1
                for i in range(l,k):
                    f = s*e[i]
                    e[i]=x*e[i]
                    if (abs(f) <= eps):
                        z = float(q[k]) 
                        if (l == k):    #masuk ke konvergen
                            if (z < 0):
                                q[k] = -z
                                for j in range(1, n):
                                    v[j,k] = -v[j,k]
                        
                        x = q[l]
                        y = q[k-1]
                        g = e[k-1]
                        h = e[k]
                        f = ((y-x)*(y+z)*(g-h)*(g+h))/(2*h*y)
                        g = math.sqrt(f*f+1)
                        if(f < 0):
                            f = ((x-z)*(x+z)+h*(y/(f-g)-h))/x
                        else:
                            f = ((x-z)*(x+z)+h*(y/(f+g)-h))/x
                        
                        c = s = 1
                        for i in range(l+1,k):
                            g = e[i]
                            y = q[i]
                            h = s*g
                            g = c*g
                            e[i-1] = z = math.sqrt(f*f+h*h)
                            c = f/z
                            s = h/z
                            f = x*c+g*s
                            g = -x*s+g*c
                            h = y*s
                            y = y*c
                            for j in range(n+1):
                                x = v[j,i-1]
                                z = v[j,i]
                                v[j,i-1] = x*c+z*s
                                v[j,i] = -x*s+z*c
                            q[i-1] = z = math.sqrt(f*f+h*h)
                            c = f/z
                            s = h/z
                            f = c*g+s*y
                            x = -s*g+c*y
                            for j in range(m+1):
                                y = u[j, i-1]
                                z = u[j,i]
                                u[j,i-1] = y*c+z*s
                                u[j,i] = -y*s+z*c
                        e[l] = 0
                        e[k] = f
                        q[k] = x
                    
                    g = q[i]
                    h = q[i] = math.sqrt(f*f+g*g)
                    c = g/h
                    s = -f/h
                    for j in range(m+1):
                        y = u[j,l1]
                        z = u[j,i]
                        u[j,l1] = y*c+z*s
                        u[j,i] = -y*s+z*c



