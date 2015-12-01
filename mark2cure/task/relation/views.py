from django.shortcuts import get_object_or_404
from django.http import HttpResponse

from django.views.decorators.http import require_http_methods
from django.contrib.auth.decorators import login_required
from django.template.response import TemplateResponse
from django.contrib import messages

from ...common.formatter import bioc_as_json
from ...common.bioc import BioCReader
from ...common.models import Group
from ...document.models import Document, Section, Pubtator

from .models import Answer, Relation, Concept
from .forms import AnswerForm


@login_required
def home(request):
    """Find documents to display to user. I was not too sure how to approach this part.For
    now, I just provide them some documents to view. This needs to be improved.
    TODO
    Just get group 4 for now. TODO upon integration, make this better (Max, ideas?)
    """
    docs_with_unanswered_relations_for_user = []

    group = Group.objects.get(pk=5)
    for document in Document.objects.filter(task__group=group).all():
        relations = Relation.objects.filter(document=document)

        for relation in relations:
            if not Answer.objects.filter(username=request.user).filter(relation=relation.pk).exists():
                if document not in docs_with_unanswered_relations_for_user:
                    docs_with_unanswered_relations_for_user.append(document)
                break

    queryset_documents = docs_with_unanswered_relations_for_user[:100]

    msg = '<p class="lead text-center">Select one of the following papers to help us further define chemical/disease relations.</p>'
    messages.info(request, msg, extra_tags='safe alert-success')

    ctx = {
        'documents': queryset_documents,
    }
    return TemplateResponse(request, 'relation/home.jade', ctx)


@login_required
def relation_api(request, document_pk):
    """Only make api for relations that do not have an answer! TODO is this okay? This can
    not be changed for now because this is why data[0] inside relation.js works. Need to
    know where the next concept is that doesn't have a user's answer already. This is
    why I do document.unanswered_relation_list (to make API for only
    unanswered_relations).
    """
    from django.http import JsonResponse
    # the document instance
    current_document = get_object_or_404(Document, pk=document_pk)

    unanswered_relations_for_user = current_document.unanswered_relation_list(request)
    # ^^^^^^ no dry TODO Max, why isn't this dry?

    relation_list = []
    pub_list = []
    for relation in unanswered_relations_for_user:
        relation_dict = {}
        concept1 = Concept.objects.get(document=current_document, uid=relation.concept1_id)
        concept2 = Concept.objects.get(document=current_document, uid=relation.concept2_id)
        relation_type = concept1.stype + "_" + concept2.stype

        def get_pub_pk(concept_stype, uid):
            global pub_list

            if concept_stype == "g":
                pub = Pubtator.objects.get(document=current_document, kind="GNormPlus")
            elif concept_stype == "d":
                pub = Pubtator.objects.get(document=current_document, kind="DNorm")
            else:
                pub = Pubtator.objects.get(document=current_document, kind="tmChem")

            if concept_stype == "g":
                pub_type = "NCBI Gene"
            elif concept_stype == "d":
                pub_type = "MEDIC"
            elif concept_stype == "c":
                pub_type = "MESH"

            pubtator = pub

            if pubtator.valid():
                pub_dict = {}
                r = BioCReader(source=pubtator.content)
                r.read()
                for doc_idx, document in enumerate(r.collection.documents):
                    for passage_idx, passage in enumerate(document.passages):
                        for annotation in r.collection.documents[doc_idx].passages[passage_idx].annotations:
                            try:
                                text = annotation.text
                                uid = annotation.infons[pub_type]
                                location = str(annotation.locations[0])
                            except:
                                continue

                            if uid is not None:
                                if uid in pub_dict:
                                    if location not in pub_dict[uid]['location']:
                                        pub_dict[uid]['location'].append(location)
                                    if text not in pub_dict[uid]['text']:
                                        pub_dict[uid]['text'].append(text)
                                else:
                                    pub_dict[uid] = {
                                        'text': [text],
                                        'stype': concept_stype,
                                        'location': [location],
                                        'uid': uid
                                    }

            return pub.pk, pub_dict

        if (relation_type == 'g_d'):
            file_key = 'gene_disease_relation_menu'

        if (relation_type == 'c_d'):
            file_key = 'chemical_disease_relation_menu'

        if (relation_type == 'g_c'):
            file_key = 'gene_chemical_relation_menu'

        # when I passed this information to html, I could not retain the object, (error "is not JSON serializable")
        # so I put them in different properties
        pub1_pk, pub1_dict = get_pub_pk(str(concept1.stype), concept1.uid)
        pub2_pk, pub2_dict = get_pub_pk(str(concept2.stype), concept2.uid)

        if pub1_pk not in pub_list:
            pub_list.append(pub1_pk)
        if pub2_pk not in pub_list:
            pub_list.append(pub2_pk)

        relation_dict['c1_id'] = str(concept1.id)
        relation_dict['c1_text'] = str(concept1.text)
        relation_dict['c1_stype'] = str(concept1.stype)
        relation_dict['c1_uid'] = str(concept1.uid)
        relation_dict['c1_locations'] = pub1_dict[concept1.uid]['location']
        relation_dict['c1_pub_pk'] = str(pub1_pk)
        relation_dict['c2_id'] = str(concept2.id)
        relation_dict['c2_text'] = str(concept2.text)
        relation_dict['c2_stype'] = str(concept2.stype)
        relation_dict['c2_uid'] = str(concept2.uid)
        relation_dict['c2_locations'] = pub2_dict[concept2.uid]['location']
        relation_dict['c2_pub_pk'] = str(pub2_pk)
        relation_dict['relation_type'] = file_key
        relation_dict['pk'] = relation.pk
        relation_dict['pub_list'] = pub_list

        relation_list.append(relation_dict)

    print relation_list
    return JsonResponse(relation_list, safe=False)


@login_required
def relation(request, document_pk):
    """Main page for users to find the relationship between two concepts."""
    document = get_object_or_404(Document, pk=document_pk)

    relations = document.unanswered_relation_list(request)
    relation = relations[0]  # was used...TODO check on this

    ctx = {'sections': Section.objects.filter(document=document),
           'relation': relation,
           'document_pk': document_pk,
           }
    return TemplateResponse(request, 'relation/relation.jade', ctx)


# Pass in relation similar to above method
def results(request, relation_id):
    # (TODO) Why is relation being passed in but never used? -Max
    relation = get_object_or_404(Relation, pk=relation_id)

    relation_type = request.POST['relation_type']
    form = AnswerForm(request.POST or None)

    # TODO: Need help with this. Can I submit *** ALL*** answers via ajax?  Do I even need
    # the final submit/form?  Seems unnecessary.
    if form.is_valid():
        print "VALID FORM"
        form.save()
    relation_type = "testing views relation_type"
    if request.method == 'POST':
        ctx = {
            'relation_type':relation_type,
        }
        return TemplateResponse(request, 'relation/results.jade', ctx)
    # return HttpResponseRedirect(reverse("relation:home"))


@login_required
@require_http_methods(['POST'])
def create_post(request):
    form = AnswerForm(request.POST or None)
    # print form
    if form.is_valid():
        form.save()

    return HttpResponse(200)


@login_required
@require_http_methods(['POST'])
def jen_bioc(request):
    # (TODO) Why are these assigned but never used? -Max
    relation_list = request.POST['relation_list']
    writer_json = bioc_as_json(relation_list)
    return HttpResponse(200)
