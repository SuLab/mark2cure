from mark2cure.common.bioc import BioCWriter, BioCCollection
import datetime


def bioc_writer(request):
    writer = BioCWriter()
    writer.collection = BioCCollection()
    writer.collection.date = datetime.date.today().strftime("%Y%m%d")
    writer.collection.source = 'Mark2Cure API: {relative_url}'.format(
        relative_url=request.META.get('PATH_INFO', ''))
    return writer
