from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from django.conf import settings

from mark2cure.common.models import Task, Group, DocumentQuestRelationship, UserQuestRelationship, SkillBadge
from mark2cure.document.models import Document, View, Annotation
from mark2cure.document.tasks import get_pubmed_document
from mark2cure.userprofile.models import UserProfile

from optparse import make_option
from brabeion import badges
import requests
import random
import csv

            user, created = User.objects.get_or_create(username='Doc_G-man')
            if created:
                UserProfile.objects.create(user=user)

                gold_profile = user.userprofile
                gold_profile.quote = "Be the change you wish to see in the world."
                gold_profile.save()

                user.set_password('')
                user.save()

                # Assign the GM User the hightest skill possible
                for index, level in enumerate(SkillBadge.levels):
                    badges.possibly_award_badge("skill_awarded", user=user, level=index)

            # Clean out all the old annotations just b/c we don't know what they were off on / need to be changed
            for doc in Document.objects.all():
                views = View.objects.filter(section__document=doc, user=user)
                for view in views:
                    Annotation.objects.filter(view=view).delete()

            gm_documents_ids.extend(mixed_gm)
            for dataset in datasets:
                with open('assets/datasets/{dataset}_annos.txt'.format(dataset=dataset), 'r') as f:
                    reader = csv.reader(f, delimiter='\t')
                    next(reader, None)  # skip the headers
                    for doc_id, doc_field, ann_type, text, start, stop in reader:
                        if int(doc_id) in gm_documents_ids:

                            doc = Document.objects.get(document_id=doc_id)
                            dqr = DocumentQuestRelationship.objects.filter(document=doc, task__experiment=settings.EXPERIMENT).first()
                            print "Document ID:", doc.pk, "Quest Relationship:", dqr.pk
                            return

                            # If the GM Anns are for a document that isn't included in any quests,
                            # don't bother
                            if dqr:
                                # Be (though uncomplete) associted with each first quest to link Views
                                gm_quest_rel, gm_quest_rel_created = UserQuestRelationship.objects.get_or_create(
                                    task=dqr.task,
                                    user=user)

                                for section in doc.section_set.all():
                                    # Make sure the annotion is for the title or abstract (our supported section types)
                                    if section.kind == doc_field[0]:
                                        view, created = View.objects.get_or_create(
                                            section=section,
                                            user=user,
                                            completed=True)

                                        gm_quest_rel.views.add(view)
                                        Annotation.objects.create(
                                            view=view,
                                            text=text,
                                            start=start,
                                            type=ann_type,
                                            kind='e')



