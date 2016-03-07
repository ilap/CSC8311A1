from django.shortcuts import render
from django.views.generic import View

from Bio import SeqIO

from .forms import *
from .models import *
from .utils import GuideRNABioUtils


class GuideRNAView(View):
    form_class = SpeciesForm
    template_new = 'grna/grna_new.html'
    template_results = 'grna/grna_results.html'

    def __init__(self):
        super(View, self).__init__()
        self.grna_utils = GuideRNABioUtils()

    def get(self, request):
        print "GET request"
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
            print "DO stuffs"
            for i in request.POST:
                print "POOOOST:", i, request.POST[i]

            species_id = request.POST['species']
            target = request.POST['target']
            pam_id = request.POST['pam']
            upstream = request.POST['upstream']
            downstream = request.POST['downstream']
            is_nickase = request.POST['is_nickase']

            # Create target model from the seqeuence
            target_model = Target()
            target_model.init_target(target, upstream, downstream)
            species_model = Species.objects.get(pk=species_id)

            grna_utils = GuideRNABioUtils()


            print "MODEL:", species_model.fasta_file

            content = grna_utils.run_blat(species_model.get_fasta_file(),
                                          target_model.get_fasta_file())

            print "CONTENT: ", content
            return render(request, self.template_results, {'content':
                                                               content, })


def grna_results(request):
    # if request.method == "POST":
    print "GRNA RESULTS...."
    for i in request.POST:
        print "POOOOST:", i, request.POST[i]

    # import time
    # time.sleep (60)
    return render(request, 'grna/grna_results.html', {})
