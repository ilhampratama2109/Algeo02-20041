from django.core.files.storage import FileSystemStorage
from .settings import MEDIA_URL
import math
from PIL import Image
import numpy as np
import io

def SVD(a) :

# KAMUS LOKAL
# i, j, k, l , l1: integer
# c, f, g, h, s, x, y, z: float
# q, e: array

# ALGORITMA

    eps = 1.e-15
    tol = 1.e-64 / eps
    assert 1.0 + eps > 1.0
    assert tol > 0.0
    m = len(a)
    n = len(a[0])

    # inisialisasi matriks 

    q = [0.0 for i in range (n)]
    u = [[0.0 for j in range (n)] for i in range (m)]
    v = [[0.0 for j in range (n)] for i in range (n)]
    e = [0.0 for i in range (n)]
    for i in range (m) :
        for j in range (n) :
            u[i][j] = a[i][j]

    # Householder's reduction to bidiagonal form

    g = x = 0.0
    for i in range (n) :
        e[i] = g
        s = 0.0
        l = i+1
        for j in range (i,m) :
            s = s+(u[j][i])**2
        if s <= tol :
            g = 0.0
        else :
            f = u[i][i]
            if (f < 0.0) :
                g = math.sqrt(s)
            else :
                g = -math.sqrt(s)
            h = f*g-s
            u[i][i] = f-g
            for j in range (l,n) :
                s = 0.0
                for k in range (i,m) :
                    s = s+u[k][i]*u[k][j]
                f = s/h
                for k in range (i,m) :
                    u[k][j] = u[k][j] + f*u[k][i]
        q[i] = g
        s = 0.0
        for j in range (l,n) :
            s = s+(u[i][j])**2
        if (s <= tol) :
            g = 0.0
        else :
            f = u[i][i+1]
            if (f < 0.0) :
                g = math.sqrt(s)
            else :
                g = -math.sqrt(s)
            h = f*g-s
            u[i][i+1] = f-g
            for j in range (l,n) :
                e[j] = u[i][j]/h
            for j in range (l,m) :
                s = 0.0
                for k in range (l,n) :
                    s = s+u[j][k]*u[i][k]
                for k in range (l,n) :
                    u[j][k] = u[j][k] + s*e[k]
        y = abs(q[i])+abs(e[i])
        if (y > x) :
            x = y

# COMMENT 3
# accumulation of right-hand transformations

    for i in range(n-1, -1, -1):
        if (g != 0.0):
            h = u[i][i+1] * g
            for j in range(l, n):
                v[j][i] = u[i][j]/h
            for j in range(l, n):
                s = 0.0
                for k in range(l, n):
                    s = s + u[i][k] * v[k][j]
                for k in range(l, n):
                    v[k][j] = v[k][j] + s * v[k][i]
        for j in range(l, n):
            v[j][i] = 0.0
            v[i][j] = v[j][i] 
        v[i][i] = 1.0
        g = e[i]
        l = i    

# COMMENT 4
# accumulation of left-hand transformations

    for i in range(n-1, -1, -1):
        l = i + 1
        g = q[i]
        for j in range(l, n):
            u[i][j] = 0.0
        if (g != 0.0):
            h = u[i][i] * g 
            for j in range(l, n):
                s = 0.0
                for k in range(l, m):
                    s = s + u[k][i] * u[k][j]
                f = s/h
                for k in range(i, m):
                    u[k][j] = u[k][j] + f * u[k][i]
            for j in range(i, m):
                u[j][i] = u[j][i]/g
        else :
            for j in range(i, m):
                u[j][i] = 0.0
        u[i][i] = u[i][i] + 1.0

# COMMENT 5
# diagonalization of tile bidiagonal form

    eps = eps*x
    for k in range(n-1, -1, -1):
        while True:
            # test f splitting
            for l in range(k,-1,-1):
                test_f_convergence = False
                if (abs(e[l]) <= eps):
                    # goto test f convergence
                    test_f_convergence = True
                    break  
                if (abs(q[l-1]) <= eps):
                    # goto cancellation
                    break 
            if not test_f_convergence:
                #cancellation of e[l] if l>0
                c = 0.0
                s = 1.0
                l1 = l-1
                for i in range(l,k+1):
                    f = s*e[i]
                    e[i] = c*e[i]
                    if abs(f) <= eps:
                        #goto test f convergence
                        break
                    g = q[i]
                    h = modSqrt(f,g)
                    q[i] = h
                    c = g/h
                    s = -f/h
                    for j in range(m):
                        y = u[j][l1]
                        z = u[j][i]
                        u[j][l1] = y*c+z*s
                        u[j][i] = -y*s+z*c
            # test f convergence
            z = q[k]
            if l == k:
                # convergence
                if z<0.0:
                    #q[k] is made non-negative
                    q[k] = -z
                    for j in range(n):
                        v[j][k] = -v[j][k]
                break  
            # shift from bottom 2x2 minor
            x = q[l]
            y = q[k-1]
            g = e[k-1]
            h = e[k]
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
            g = modSqrt(f,1.0)
            if f < 0:
                f = ((x-z)*(x+z)+h*(y/(f-g)-h))/x
            else:
                f = ((x-z)*(x+z)+h*(y/(f+g)-h))/x
            # next QR transformation
            c = 1.0
            s = 1.0
            for i in range(l+1,k+1):
                g = e[i]
                y = q[i]
                h = s*g
                g = c*g
                z = modSqrt(f,h)
                e[i-1] = z
                c = f/z
                s = h/z
                f = x*c+g*s
                g = -x*s+g*c
                h = y*s
                y = y*c
                for j in range(n):
                    x = v[j][i-1]
                    z = v[j][i]
                    v[j][i-1] = x*c+z*s
                    v[j][i] = -x*s+z*c
                z = modSqrt(f,h)
                q[i-1] = z
                c = f/z
                s = h/z
                f = c*g+s*y
                x = -s*g+c*y
                for j in range(m):
                    y = u[j][i-1]
                    z = u[j][i]
                    u[j][i-1] = y*c+z*s
                    u[j][i] = -y*s+z*c
            e[l] = 0.0
            e[k] = f
            q[k] = x
            # goto test f splitting

    vt = transposeMatrix(v)
    return u,q,vt

def transposeMatrix(M):
    m = len(M)
    n = len(M[0])
    TM = []
    for i in range(n): TM.append([0.0]*m) # inisiliasi matriks
    for i in range(m):
        for j in range(n):
            TM[j][i] = M[i][j]
    return TM

def modSqrt(a, b):
    absa = abs(a)
    absb = abs(b)
    if absa > absb: 
        return absa*math.sqrt(1.0+(absb/absa)**2)
    else:
        if absb == 0.0: 
            return 0.0
        else: 
            return absb*math.sqrt(1.0+(absa/absb)**2)

def compress(image_url, k):
    k = int(k)
    # image processing
    img = Image.open("C:/Tubes Algeo/Tubes 2/Algeo02-20041/src" + image_url)
    img_format = img.format
    print(format)
    image = np.array(img)
    image = image/255
    row,col,_ = image.shape

    print("pixels: ",row, " ", col)

    image_red = image[:,:,0]
    image_green = image[:,:,1]
    image_blue = image[:,:,2]

    landscape = False
    if col > row:
        landscape = True
        image_red = np.transpose(image_red)
        image_green = np.transpose(image_green)
        image_blue =  np.transpose(image_blue)

    u_r,q_r,v_r = SVD(image_red)
    u_g,q_g,v_g = SVD(image_green)
    u_b,q_b,v_b = SVD(image_blue)

    urk = np.array(u_r)[:,0:k]
    vrk = np.array(v_r)[0:k,:]
    ugk = np.array(u_g)[:,0:k]
    vgk = np.array(v_g)[0:k,:]
    ubk = np.array(u_b)[:,0:k]
    vbk = np.array(v_b)[0:k,:]
    qrk = np.array(q_r)[0:k]
    qgk = np.array(q_g)[0:k]
    qbk = np.array(q_b)[0:k]

    image_red_approx = np.dot(urk,np.dot(np.diag(qrk),vrk))
    image_green_approx = np.dot(ugk,np.dot(np.diag(qgk),vgk))
    image_blue_approx = np.dot(ubk,np.dot(np.diag(qbk),vbk))

    if landscape:
        image_red_approx = np.transpose(image_red_approx)
        image_green_approx = np.transpose(image_green_approx)
        image_blue_approx = np.transpose(image_blue_approx)

    image_recons = np.zeros((row,col,3))
    image_recons[:,:,0] = image_red_approx
    image_recons[:,:,1] = image_green_approx
    image_recons[:,:,2] = image_blue_approx

    image_recons[image_recons < 0] = 0
    image_recons[image_recons > 1] = 1

    im = Image.fromarray(np.uint8(image_recons*255))
    b = io.BytesIO() # "convert" image object to image file
    im.save(b, img_format)
    fs = FileSystemStorage()
    compressed_img = fs.save("hasil." + img_format, b)
    b.seek(0)
    #im.save("hasil.jpg")
    print("Selesai")
    return fs.url(compressed_img)
    
