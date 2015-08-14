


from **TASKS**, which is working
(minus lingpipe and sentence splitter). **models** is not working
entirely yet.

using python manage.py shell_plus, execute the following commands
to work with the "relation module features"

from mark2cure.relation.tasks import parse_input

import os

papers = parse_input(os.getcwd(),"CDR_small.txt")

In [6]: print papers

output:
[<PaperTask>: PMID 227508. 20 annotations, 0 gold relations
3 unique chemical ids, 2 unique disease ids
2 sentences, <PaperTask>: PMID 354896. 6 annotations, 0 gold relations
1 unique chemical ids, 3 unique disease ids
2 sentences, <PaperTask>: PMID 435349. 14 annotations, 0 gold relations
1 unique chemical ids, 2 unique disease ids
2 sentences]


papers[0].get_unique_concepts()
output is **unique** chemical and disease IDs. There are many repetitive annotations, so we want to reduce the amount of work to make things more efficient for volunteers.
({'D003000', 'D008750', 'D009270'}, {'D006973', 'D007022'})



NEW, working html is verify_relation.html
