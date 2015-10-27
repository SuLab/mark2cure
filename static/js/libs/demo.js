var data = {
  "chemical_disease_relation_menu": [{
    "id": "1",
    "text": "(chemical) relates to disease",
    "desc": "The chemical has some type of relation to the disease.",
    "example": "Evidence indicates that Amblyomin-X could be a promising candidate for cancer therapy.",
    "children": [{
        "id": "1-2",
        "text": "Exacerbates",
        "desc": "The chemical worsens or exacerbates the disease in any manner.",
        "example": "St. John's Wort has been known to worsen symptoms of Parkinson's disease.",
      }, {
        "id": "1-1",
        "text": "Treats",
        "desc": "The chemical improves or treats the disease in any manner.",
        "example": "Treatment with gingerol may provide a novel therapeutic strategy for prion-mediated neurotoxicity",
      }, {
        "id": "1-3",
        "text": "Increases risk of disease",
        "desc": "Administration of the drug increases chances of getting the disease. Considered an 'adverse side effect' of the drug.",
        "example": "Parkinson's disease risk factors include ageing, exposure to toxins and genetic defects.",
      }, {
        "id": "1-4",
        "text": "May cause",
        "desc": "Administration of the chemical may cause or lead to the disease.",
        "example": "43% of mice given the chemical NVP-BHG712 developed pulmonary metastases.",
      }, {
        "id": "1-5",
        "text": "Prevents",
        "desc": "Administration of the chemical prevents the disease",
        "example": "Malarone is taken daily by travelers wanting to prevent malaria.",
      }]
    }, {
    "id": "2",
    "text": "has no relation to",
    "desc": "There is no relation between the chemical and disease.",
    "example": "To study soft tissue sarcoma, an experimental group of mice was given the chemical NVP-BHG712.",
  }],

  "gene_disease_relation_menu": [{
    "id": "1",
    "text": "(gene) relates to the disease",
    "desc": "The gene is related to the disease in some manner.",
    "example": "The gene responsible for Triple A syndrome, AAAS, has recently been identified.",
    "children": [{
      "id": "1-1",
      "text": "Altered expression is associated with",
      "desc": "Gene expression is the process by which information from a gene is used in the synthesis of a functional gene product (mRNA, protein, or microRNAs). When gene expression is altered, gene expression is either increased or decreased.",
      "example": "Several studies revealed significantly higher EPHB4 expression in malignancies such as prostate cancer.",
      }, {
      "id": "1-3",
      "text": "Mutation is associated with",
      "desc": "Gene mutations or aberrations are permanent changes in the gene DNA sequence, differing from the sequence found in most people. Mutations range in size; they can affect anywhere from a single DNA base pair to a large piece of a chromosome that includes multiple genes.",
      "example": "Mutations in the COL5A or COL3A genes are only a few of the genetic causes of Ehlers-Danlos syndrome.",
      }, {
        "id": "1-4",
        "text": "And post-translational modifications are associated with",
        "desc": "Translation is protein synthesis. When changes are made to a protein during or after translation, it is considered a post-translational modification.",
        "example": "A small interfering RNA causes knockdown of ATP2C1 expression, resulting in defects in both post-translational processing of wild-type thyroglobulin",
    }]}, {
    "id": "2",
    "text": "Has no relation to",
    "desc": "The gene does not relate to the disease.",
    "example": "The precise role of Ngly1 in the ERAD process remains unclear in mammals.",
  }],

  "gene_chemical_relation_menu": [{
    "id": "1",
    "text": "(gene) relates to the chemical",

    "children": [{
      "id": "1-1",
      "text": "amount is changed by",

      "children": [{
        "id": "1-1-1",
        "text": "expression is decreased by"
      }, {
        "id": "1-1-2",
        "text": "expression is increased by"
      }]

    }, {
      "id": "1-2",
      "text": "activity is changed by",

      "children": [{
        "id": "1-2-1",
        "text": "activity is inhibited by"
      }, {
        "id": "1-2-2",
        "text": "activity is activated by"
      }]
    }]

    }, {
    "id": "2",
    "text": "has no relation to"
    }
  ]

};


Tree.addInitializer(function(options) {

  /* data fetches the original document and annotations
  $.getJSON('http://localhost:12345/data/relationships.json', function( data ) {
  */
    Tree.addRegions({'start': '#tree-insert'});
    console.log("testing inside demo.js");
    /* When the app is first loaded */
    var coll = new RelationList(data[file_key]);

    Tree['start'].show( new RelationCompositeView({
      collection: coll,
      concepts: concepts,
      choice: false
    }));

    Backbone.Radio.DEBUG = true;
    Tree['convoChannel'] = Backbone.Radio.channel('convo');

    /* When an item is selected */
    Tree['convoChannel'].on('click', function(obj) {
      Tree['start'].show( new RelationCompositeView({
        collection: obj['collection'],
        concepts: concepts,
        choice: obj['choice']
      }));
      console.log(obj['choice'].id);
    });

    /* When the back toggle is selected */
    Tree['convoChannel'].on('back', function(opts) {
      var collection;

      var parentRel = opts['choice'].get('parentRelation');
      if(parentRel) { collection = parentRel.get('children'); }

      if(opts['collection']) { collection = opts['collection']; }

      /* Backup: Go to the top of the stack */
      collection = coll;

      /* Call the View Redraw */
      Tree['start'].show( new RelationCompositeView({
        collection: coll,
        concepts: concepts,
        choice: false
      }));

    });

  /*
  });
  */
});

Tree.start();

/* Other stuff for Rel demo */
$('.entity').on('mouseover', function(evt) {
  $(this).find('.text').hide();
  $(this).find('.message').show();
});

$('.entity').on('mouseout', function(evt) {
  $(this).find('.text').show();
  $(this).find('.message').hide();
});
