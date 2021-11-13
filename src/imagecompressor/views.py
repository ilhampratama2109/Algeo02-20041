from django.shortcuts import render
from django.core.files.storage import FileSystemStorage
from .forms import ImageForm
from .svd import compress

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
        print(uploaded_image.name)
        fs = FileSystemStorage()
        name = fs.save(uploaded_image.name, uploaded_image)
        print(fs.url(name))
        context['img_url'] = fs.url(name)
        context['rate'] = request.POST['rate']
        context['img_compressed_url'], context['runtime'], context['percentage']  = compress(context['img_url'], request.POST['rate'])

    return render(request, "index.html", context)