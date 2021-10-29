from django.shortcuts import render

def index(request):
    context = {
        'heading' : 'Image Compressor',
        'title' : 'Image Compressor',
    }
    return render(request, "index.html", context)