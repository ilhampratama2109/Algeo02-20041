from django import forms 
from django.core.validators import MaxValueValidator, MinValueValidator

class ImageForm(forms.Form):
    image = forms.ImageField(
        label='Choose your image',
        widget = forms.FileInput(
            attrs={
                'class':'form-control'
            }
        )
    )
    rate = forms.IntegerField(
        min_value=1,
        label='Image compression rate', 
        widget = forms.NumberInput(
            attrs={
                'class':'form-control',
                'placeholder':'Input rate (%)',
            }
        )    
    )
