from django.shortcuts import render

# Create your views here.

def index(request):
    context = {
        'heading' : 'Tentang Kami',
        'title' : 'Tentang Kami',
        'nav' : [
            ['/', 'Beranda'],
        ],
    }
    return render(request, 'about/index.html', context)
