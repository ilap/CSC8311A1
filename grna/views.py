from django.shortcuts import render

def grna_new (request):
    return render(request, 'grna/grna_new.html', {})

def grna_results (request):
    if request.method == "POST":

        for i in request.POST:
            print "POOOOST:", i, request.POST[i]

    #import time
    #time.sleep (60)
    return render(request, 'grna/grna_results.html', {})


