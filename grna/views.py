from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse

from .forms import *
from .models import *


def grna_new(request):

    if request.POST.get('upload'):
        print "Upload"
        form = SpeciesForm(request.POST, request.FILES)
        if form.is_valid():
            # file is saved
            species = Species(fasta_file=request.FILES['fasta_file'])
            species.save()

            # Redirect to the document list after POST
            return HttpResponseRedirect(reverse('myapp.views.list'))
        return render(request, 'grna/grna_new.html', {})

    elif request.POST.get('search'):
        print "DO stuffs"
        return render(request, 'grna/grna_results.html', {})

    else:
        form = SpeciesForm()

        species = Species.objects.all()
        nucleases = Nuclease.objects.all()
        pams = PAM.objects.all()

        print species
        print nucleases
        print pams

        return render(request, 'grna/grna_new.html', {'form': form,
                                                  'species': species,
                                                  'nucleases': nucleases,
                                                  'pams': pams, })


def grna_results(request):
    if request.method == "POST":

        for i in request.POST:
            print "POOOOST:", i, request.POST[i]

    # import time
    # time.sleep (60)
    return render(request, 'grna/grna_results.html', {})
