from django import template

# Filtering
register = template.Library()


@register.filter
def in_hits(grnas, hit):
    return grnas.filter(target_hit=hit)


@register.filter
def percentage(value):
    return '{0:.2%}'.format(value)
