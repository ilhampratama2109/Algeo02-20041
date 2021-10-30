from django.shortcuts import render

def index(request):
    context = {
        'heading' : 'Image Compressor',
        'title' : 'Image Compressor',
        'nav' : [
            ['/', 'Beranda'],
            ['/about', 'Tentang Kami'],
        ]
    }
    return render(request, "index.html", context)