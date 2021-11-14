from django.core.files.storage import FileSystemStorage
from PIL import Image
from time import time
import numpy as np
import io
import math

def SVD(a) :

# KAMUS LOKAL
# i, j, k, l, g, l1: integer
# c, f, g, h, s, x, y, z: float
# q, e: array

# ALGORITMA

    eps = 1.e-15 # konstanta yang digunakan dalam tes konvergensi
    tol = 1.e-64 / eps # konstanta yang bergantung kepada komputer/mesin yg digunakan
    assert 1 + eps > 1
    assert tol > 0
    m = len(a) # mendapatkan jumlah kolom
    n = len(a[0]) # mendapatkan jumlah baris

    # inisialisasi matriks 

    q = [0 for i in range (n)]
    u = [[0 for j in range (n)] for i in range (m)]
    v = [[0 for j in range (n)] for i in range (n)]
    e = [0 for i in range (n)]
    for i in range (m) :
        for j in range (n) :
            u[i][j] = a[i][j]

    # reduksi ke dalam bentuk bidiagonal

    g = x = 0
    for i in range (n) :
        e[i] = g
        s = 0
        l = i+1
        for j in range (i,m) :
            s = s+(u[j][i])**2
        if s <= tol :
            g = 0
        else :
            f = u[i][i]
            if (f < 0) :
                g = math.sqrt(s)
            else :
                g = -math.sqrt(s)
            h = f*g-s
            u[i][i] = f-g
            for j in range (l,n) :
                s = 0
                for k in range (i,m) :
                    s = s+u[k][i]*u[k][j]
                f = s/h
                for k in range (i,m) :
                    u[k][j] = u[k][j] + f*u[k][i]
        q[i] = g
        s = 0
        for j in range (l,n) :
            s = s+(u[i][j])**2
        if (s <= tol) :
            g = 0
        else :
            f = u[i][i+1]
            if (f < 0) :
                g = math.sqrt(s)
            else :
                g = -math.sqrt(s)
            h = f*g-s
            u[i][i+1] = f-g
            for j in range (l,n) :
                e[j] = u[i][j]/h
            for j in range (l,m) :
                s = 0
                for k in range (l,n) :
                    s = s+u[j][k]*u[i][k]
                for k in range (l,n) :
                    u[j][k] = u[j][k] + s*e[k]
        y = abs(q[i])+abs(e[i])
        if (y > x) :
            x = y

    # A.A(transpose)

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

    # A(transpose).A

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


    # diagonalisasi bentuk bidiagonal

    eps = eps*x
    for k in range(n-1, -1, -1):
        while True:
            # tes f splitting
            for l in range(k,-1,-1):
                test_f_convergence = False
                if (abs(e[l]) <= eps):
                    # masuk ke test f convergence
                    test_f_convergence = True
                    break  
                if (abs(q[l-1]) <= eps):
                    # masuk ke cancellation
                    break 
            if not test_f_convergence:
                # cancellation 
                c = 0
                s = 1
                l1 = l-1
                for i in range(l,k+1):
                    f = s*e[i]
                    e[i] = c*e[i]
                    if abs(f) <= eps:
                        # masuk ke test f convergence
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
            if (l == k):
                # convergence
                if (z < 0):
                    # q[k] dibuat tidak negatif
                    q[k] = -z
                    for j in range(n):
                        v[j][k] = -v[j][k]
                break  
            x = q[l]
            y = q[k-1]
            g = e[k-1]
            h = e[k]
            f = ((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y)
            g = modSqrt(f,1)
            if (f < 0):
                f = ((x-z)*(x+z)+h*(y/(f-g)-h))/x
            else:
                f = ((x-z)*(x+z)+h*(y/(f+g)-h))/x
            # transformasi QR selanjutnya
            c = 1
            s = 1
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
            e[l] = 0
            e[k] = f
            q[k] = x
            # masuk ke tes f splitting

    vt = transposeMatrix(v)
    return u,q,vt

def transposeMatrix(M): # fungsi transpose matriks 
    # KAMUS LOKAL
    # m, n: integer
    # TM: matriks

    # ALGORITMA
    m = len(M)
    n = len(M[0])
    TM = []
    for i in range(n): TM.append([0]*m) # inisiliasi matriks
    for i in range(m):
        for j in range(n):
            TM[j][i] = M[i][j]
    return TM

def modSqrt(a, b): # fungsi sqrt yang dimodifikasi
    # KAMUS LOKAL
    # x, y : integer

    # ALGORITMA
    x = abs(a)
    y = abs(b)
    if (x > y): 
        return x*math.sqrt(1+(y/x)**2)
    else:
        if (y == 0): 
            return 0
        else: 
            return y*math.sqrt(1+(x/y)**2)

def compress(image_url, percent): # fungsi kompresi gambar
    start = time()

    img = Image.open("C:/Tubes Algeo/Tubes 2/Algeo02-20041/src" + image_url) # disesuaikan dengan directory projek
    img_format = img.format
    print(format)
    image = np.array(img)
    image = image / 255
    row,col,_ = image.shape
    k = round(int(percent) * row * col / (100 * (row + col + 1)))

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
    b = io.BytesIO() # "convert" objek Image ke bentuk file gambar
    im.save(b, img_format)
    fs = FileSystemStorage()
    compressed_img = fs.save("hasil." + img_format, b)
    b.seek(0)

    end = time()
    run_time = "{0:.2f}".format(end - start)
    print("Selesai")
    return fs.url(compressed_img), run_time, k
