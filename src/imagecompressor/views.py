from django.shortcuts import render
from django.core.files.storage import FileSystemStorage
from .forms import ImageForm

def index(request):
    image_form = ImageForm()
    context = {
        'heading' : 'Image Compressor',
        'title' : 'Image Compressor',
        'nav' : [
            ['/', 'Beranda'],
            ['/about', 'Tentang Kami'],
        ],
        'data_form' : image_form,
    }

    if request.method == 'POST':
        uploaded_image = request.FILES['image']
        fs = FileSystemStorage()
        name = fs.save(uploaded_image.name, uploaded_image)
        print(fs.url(name))
        context['img_url'] = fs.url(name)
    return render(request, "index.html", context)