from django.shortcuts import render

# Create your views here.
from django.http import HttpResponse
from django.template import loader


from .models import GuideRNA


def index(request):
    return HttpResponse("Hello, world. You're at the polls index.")

def iindex(request):
    latest_question_list = GuideRNA.objects.order_by('-pub_date')[:5]
    template = loader.get_template('grna/index.html')
    context = {
        'latest_question_list': latest_question_list,
    }
    return HttpResponse(template.render(context, request))

def indea2x(request):
    latest_question_list = GuideRNA.objects.order_by('-pub_date')[:5]
    output = ', '.join([q.question_text for q in latest_question_list])
    return HttpResponse(output)

def detail(request, question_id):
    return HttpResponse("You're looking at question %s." % question_id)

def results(request, question_id):
    response = "You're looking at the results of question %s."
    return HttpResponse(response % question_id)

def vote(request, question_id):
    return HttpResponse("You're voting on question %s." % question_id)