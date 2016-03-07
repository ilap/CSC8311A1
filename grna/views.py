from django.shortcuts import render
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.views.generic import View

from .forms import *
from .models import *

class GuideRNAView(View):
    form_class = SpeciesForm
    template_new = 'grna/grna_new.html'
    template_results = 'grna/grna_results.html'

    def get(self, request, *args, **kwargs):
        form = self.form_class()
        species = Species.objects.all()
        nucleases = Nuclease.objects.all()
        pams = PAM.objects.all()

        print species
        print nucleases
        print pams
        return render(request, self.template_new, {'form': form,
                                              'species': species,
                                              'nucleases': nucleases,
                                              'pams': pams, })

    def post(self, request, *args, **kwargs):
        if request.POST.get('upload'):
            print "Upload"
            form = self.form_class(request.POST, request.FILES)
            if form.is_valid():
                # file is saved
                species = Species(fasta_file=request.FILES['fasta_file'])
                species.save()

                # Redirect to the document list after POST
                return HttpResponseRedirect(reverse('myapp.views.list'))
            return render(request, self.template_new, {})
        elif request.POST.get('search'):
            print "DO stuffs"
            return render(request, self.template_results, {})


def grna_results(request):
    if request.method == "POST":

        for i in request.POST:
            print "POOOOST:", i, request.POST[i]

    # import time
    # time.sleep (60)
    return render(request, 'grna/grna_results.html', {})
