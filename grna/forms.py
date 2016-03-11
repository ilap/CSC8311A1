from django import forms


# Not used currently, it will be used for uploading files.
class SpeciesForm(forms.Form):
    title = forms.CharField(max_length=50)
    fasta_file = forms.FileField(label='Select a file', help_text='fasta file')
