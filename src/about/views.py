from django.shortcuts import render

# Create your views here.

def index(request):
    context = {
        'heading' : 'About Us',
        'title' : 'About Us',
        'nav' : [
            ['/', 'Home'],
        ],
    }
    return render(request, 'about/index.html', context)
