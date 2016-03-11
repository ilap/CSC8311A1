from django.shortcuts import render
from django.views.generic import View
from django.shortcuts import render_to_response

from django.template import Context
from django.template import RequestContext

from .forms import *
from .models import *
from .utils import GuideRNAManager

class GuideRNAView(View):
    form_class = SpeciesForm
    template_new = 'grna/grna_new.html'
    template_results = 'grna/grna_results.html'

    def __init__(self):
        super(View, self).__init__()
        self.grna_utils = GuideRNAManager()


    def get(self, request):
        print "GET request"
        form = self.form_class()
        species = Species.objects.all()
        nuclease = Nuclease.objects.get(name='Cas9')
        pams = PAM.objects.filter(nuclease=nuclease)
        return render(request, self.template_new, {'form': form,
                   'species': species,
                   'nuclease': nuclease,
                   'pams': pams, })

    def post(self, request):
        print "POST request"
        if request.POST.get('upload'):
            pass
            # print "Upload"
            # form = self.form_class(request.POST, request.FILES)
            # if form.is_valid():
            #    # file is saved
            #    species = Species(fasta_file=request.FILES['fasta_file'])
            #    species.save()

            #    # Redirect to the document list after POST
            #    return HttpResponseRedirect(reverse('views.upload'))
            # return render(request, self.template_new, {})
        elif request.POST.get('search'):

            # DEBUG for i in request.POST:
            # DEBUG    print "POOOOST:", i, request.POST[i]
            # Some erro handling
            try:
                grna_utils = GuideRNAManager().initialise_run(request)
                hits = grna_utils.search_grna()
                grnas = GuideRNA.objects.all()
                species = grna_utils._species
                target = grna_utils._target
            except:
                error = True
                form = self.form_class()
                species = Species.objects.all()
                nuclease = Nuclease.objects.get(name='Cas9')
                pams = PAM.objects.filter(nuclease=nuclease)
                return render(request, self.template_new, {'form': form,
                               'species': species,
                               'nuclease': nuclease,
                               'pams': pams,
                                'error':error})




            return render(request, self.template_results,
                              {'species':species, 'hits':hits, 'grnas':grnas,
                               'target':target})

def grna_results(request):
    # import time
    # time.sleep (60)
    return render(request, 'grna/grna_results.html', {})
