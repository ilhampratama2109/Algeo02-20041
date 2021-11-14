# Algeo02-20041
Tubes 2 Algeo

Website kompresi gambar menggunakan metode Singular Value Decomposition (SVD). Algoritma perhitungan SVD yang digunakan adalah algoritma yang dikemukakan oleh G. H. Golub dan C. Reinsch.

Anggota :
1. Ilham Pratama - 13520041
2. Mohamad Daffa Argakoesoemah - 13520118
3. Haidar Ihzaulhaq - 13520150

## Tools
- Django 3.2.8
- Bootstrap 5.0
- Numpy
- Pillow

## Setup
Langkah setup dibawah dilakukan pada sistem operasi Windows. Sistem operasi lainnya silakan menyesuaikan.
1. Unduh dan install Python3 melalui https://www.python.org/downloads/
2. Masuk ke cmd. Setup virtual environment (direkomendasikan) pada directory project dengan perintah `python -m venv namafolder`. Ganti namafolder sesuai namafolder yang diinginkan
3. Di dalam cmd, masuk ke dalam directory: namafolder/Scripts lalu ketik perintah `activate`
4. Pada cmd, jika sebelum root folder terdapat `(namafolder)` artinya Anda berhasil masuk ke dalam virtual env
5. Unduh Django versi 3.2.8 dengan perintah `pip install Django==3.2.8`
6. Unduh numpy dengan perintah `pip install numpy`
7. Unduh Pillow dengan perintah `pip install pillow`
6. Cek apakah ketiga tools di atas sudah terunduh dengan perintah `pip list`

## Menjalankan Projek pada localhost
1. Masuk ke dalam virtual env dengan mengikuti langkah 3 pada bagian "Setup" di atas
2. **Ganti path untuk membuka image** pada prosedur compress dalam file svd.py sesuai path tempat Anda menyimpan projek. Petunjuk: 

   Pada baris kode di bawah:
   ``` img = Image.open("C:/Tubes Algeo/Tubes 2/Algeo02-20041/src" + image_url) ```
   
   Ganti ``` "C:/Tubes Algeo/Tubes 2/Algeo02-20041/src" ``` sesuai path tempat menyimpan projek.

3. Jalankan command `python manage.py runserver`
4. Buka di dalam browser alamat localhost yang muncul pada cmd
