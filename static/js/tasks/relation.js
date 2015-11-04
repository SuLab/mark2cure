var data = {
  "chemical_disease_relation_menu": [{
    "id": "1",
    "text": "relates to disease",
    "desc": "The chemical has some type of relation to the disease.",
    "example": "Evidence indicates that Amblyomin-X could be a promising candidate for cancer therapy.",
    "children": [{
        "id": "1-2",
        "text": "exacerbates",
        "desc": "The chemical worsens or exacerbates the disease in any manner.",
        "example": "St. John's Wort has been known to worsen symptoms of Parkinson's disease.",
      }, {
        "id": "1-1",
        "text": "treats",
        "desc": "The chemical improves or treats the disease in any manner.",
        "example": "Treatment with gingerol may provide a novel therapeutic strategy for prion-mediated neurotoxicity",
      }, {
        "id": "1-3",
        "text": "increases risk of disease",
        "desc": "Administration of the drug increases chances of getting the disease. Considered an 'adverse side effect' of the drug.",
        "example": "Parkinson's disease risk factors include ageing, exposure to toxins and genetic defects.",
      }, {
        "id": "1-4",
        "text": "may cause",
        "desc": "Administration of the chemical may cause or lead to the disease.",
        "example": "43% of mice given the chemical NVP-BHG712 developed pulmonary metastases.",
      }, {
        "id": "1-5",
        "text": "prevents",
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
    "text": "relates to the disease",
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
    "text": "relates to the chemical",

    "children": [{
      "id": "1-1",
      "text": "amount is changed by",
    }, {
      "id": "1-2",
      "text": "activity is changed by",
    }]

    }, {
    "id": "2",
    "text": "has no relation to"
    }
  ]

};

Tree.addInitializer(function(options) {

  Tree.addRegions({'start': '#tree-insert'});
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
});


/* Other stuff for Rel demo */
$('.entity').on('mouseover', function(evt) {
  $(this).find('.text').hide();
  $(this).find('.message').show();
});

$('.entity').on('mouseout', function(evt) {
  $(this).find('.text').show();
  $(this).find('.message').hide();
});


/* The current relationship pair we are working on in the array
At document start, we should work with first item, then "ON SUBMIT" we should
move through the relation_list array and bring up new concept pairs.

If a concept was flagged as incorrect, **do not** pull up that concept again,
so iterate through the array until a new concept is found for the user.
*/
function format_text_colors(section_text, section1_length, section2_length, relationship_obj, highlight_dict, section_count) {

  function color_find(relationship_type) {
    var color;
    if (relationship_type === "g") { color = "#B1FFA8"; };
    if (relationship_type === "d") { color = "#d1f3ff"; };
    if (relationship_type === "c") { color = "#ffd1dc"; };
    return color;
  };

  str = section_text;
  console.log(highlight_dict, "highlight_dict");
  console.log(section1_length, "section1_length");
  console.log(section2_length, "section2_length");
  console.log(section_text, "section_text");
  console.log(section_text.length, ".length");

  offset = 0;

  keys = Object.keys(highlight_dict).sort(function(a, b) {return a - b;});

  // process keys in reverse order
  for (i = keys.length - 1; i > -1; i--) {
    console.log(section_count, "section_count", typeof(section_count), "section_count_type");
    key = keys[i];
    console.log(key, "key");
    key_offset = key;
    if (section_count == 0 && key > section1_length) {  // Passed the end of section 1
      console.log("continue")
      continue;
    }
    else if (section_count == 1 && key < section1_length) {  // In section 2 but key is in section 1
      console.log("continue")
      continue;
    }
    else if (section_count == 1 && key > section1_length) {  // In section 2 and need to shift left by section1_length
      key_offset = key - section1_length - 1;
    }
    else {  // Within section 1
      key_offset = key;
    };
    console.log(key_offset, "key_offset");
    try {
    console.log(highlight_dict[key][0]['span'], "highlight_dict[key]");
    if (highlight_dict[key][0]['text'] == relationship_obj.c1_text){
      concept_idx = 1;
    } else {
      concept_idx = 2;
    };
    str = str.replace(new RegExp("^(.{" + key_offset + "})(.{" + highlight_dict[key][0]['span'] + "})"),   "$1<strong style='background-color:" + color_find(relationship_obj['c'+concept_idx+'_stype']) + "'>$2</strong>");

  }
  catch(err){
    console.log(err, "LOG THE ERROR");
  };

  };
  return str;
};

var global_data;
var concept_pairs_remaining;
var concept_pairs_total;
var clean_section_text;
var section_text1_offset;
var section_text2_offset;



// Makes object containing text locations
function get_annotation_spans(relationship_obj, passage1, passage2) {
  highlight_dict_LARGE = {};
  concepts = [relationship_obj.c1_text, relationship_obj.c2_text];
  console.log(concepts, "concepts");
  console.log(passage1, "passage1");
  console.log(passage2, "passage2");

  for (i = 0; i < concepts.length; i++) {
    highlight_list_small = [];
    concept = concepts[i];
    for (ann in passage1) {
      if (passage1[ann]['text'] == concept) {
        highlight_dict_small = {};
        highlight_dict_small['index'] = passage1[ann]['location']['@offset'];
        highlight_dict_small['span'] = passage1[ann]['location']['@length'];
        highlight_dict_small['text'] = passage1[ann]['text'];
        console.log(highlight_list_small, "highlight_dict_small");
        highlight_list_small.push(highlight_dict_small);
        highlight_dict_LARGE[highlight_dict_small['index']] = highlight_list_small
      };
    };
    for (ann in passage2) {
      if (passage2[ann]['text'] == concept) {
        highlight_dict_small = {};
        highlight_dict_small['index'] = passage2[ann]['location']['@offset'];
        highlight_dict_small['span'] = passage2[ann]['location']['@length'];
        highlight_dict_small['text'] = passage2[ann]['text'];
        console.log(highlight_list_small, "highlight_dict_small");
        highlight_list_small.push(highlight_dict_small);
        highlight_dict_LARGE[highlight_dict_small['index']] = highlight_list_small
      };
    };
  };
  console.log(highlight_dict_LARGE, "higlight_dict_LARGE");
  return highlight_dict_LARGE;
};


$.getJSON('/relation/'+ document_pk +'/api/', function(data) {
  global_data = data;

  // console.log(data);
  concept_pairs_remaining = data.length;
  console.log(concept_pairs_remaining, "concept pairs remaining");
  concept_pairs_total = concept_pairs_remaining;
  var relationship_obj = data[0];


// 9886 chemicals, OCP2 and 15 N
// 9888 genes, just OCP2
  $.getJSON('/document/pubtator/specific/'+ '9886' +'.json', function( data_pubtator_tmchem ) {
    pubtator_data = data_pubtator_tmchem;

    passage1 = pubtator_data.collection.document.passage[0]['annotation'];
    passage2 = pubtator_data.collection.document.passage[1]['annotation'];
    section_text1_offset = pubtator_data.collection.document.passage[0]['text'].length
    section_text2_offset = pubtator_data.collection.document.passage[1]['text'].length


    // reformat each section and use highlight_dict to change colors in the text at the same time.
    section_count = 0;
    $('.section').each(function(el_idx, el) {
      var section_text = $(this).text();
      highlight_dict = get_annotation_spans(relationship_obj, passage1, passage2);
      section_text = format_text_colors(section_text, section_text1_offset, section_text2_offset, relationship_obj, highlight_dict, section_count);
      $(this).html(section_text);
      section_count++;
    });
  });


  // Start the tree
  file_key = relationship_obj.relation_type;
  concepts = {
    'c1': {'text': relationship_obj.c1_text, 'type': relationship_obj.c1_stype },
    'c2': {'text': relationship_obj.c2_text, 'type': relationship_obj.c2_stype }
  };
  Tree.start();


});

function yays_selection() {
  yays = ['Sweet!', 'Great!', 'Nice!', 'Awesome!', 'Nice job!', 'Excellent!', 'Wow!', 'Woohoo!', 'Hooray!', 'Yeaahh!', 'Look at you!'];
  yay_word = yays[Math.floor(Math.random()*yays.length)];
  return yay_word;
};

function carry_on_selection(concept_pairs_total, concept_pairs_remaining) {
  carry_on = ['Carry on.', 'Keep going.', "Help some more.", "We still need you.", "Don't stop now.", "", "", "" ];
  if (concept_pairs_total >= 5 && concept_pairs_remaining <= 3 )
    carry_on = ["You're a trooper.", 'Most impressive.', "We know you're dedicated.", 'Very persistent.', 'Almost there.', 'On a role.', 'You are amazing.'];
  carry_on_word = carry_on[Math.floor(Math.random()*carry_on.length)];
  return carry_on_word;
};

function html_praise_display(concept_pairs_total, concept_pairs_remaining) {
  console.log(concept_pairs_total);
  console.log(concept_pairs_remaining);
  if (concept_pairs_remaining == 1) {
    html_display = "<em>" + yays_selection() + " This is the last one in this text!</em>";
  }
  else {
    html_display = "<em>" + yays_selection() + " Only " + concept_pairs_remaining + " concept pairs left. " + carry_on_selection(concept_pairs_total, concept_pairs_remaining) +"</em>";
  }
  return html_display;
};


var counter = 1;

$('#submit_button').on('click', function(evt) {
  concept_pairs_remaining -= 1;
  if (concept_pairs_remaining == 0) {
    pass;
  }
  else {
    evt.preventDefault();
  };
  $('#praise_after_completion').hide();
  var relationship_obj = global_data[counter];
  counter += 1;

  console.log(concept_pairs_remaining);
  html_display = html_praise_display(concept_pairs_total, concept_pairs_remaining);
  $('#praise_after_completion').html( html_display ).fadeIn(900);

  section_count = 0;
  $('.section').each(function(el_idx, el) {
    var section_text = $(this).text();
    highlight_dict = get_annotation_spans(relationship_obj, passage1, passage2);
    section_text = format_text_colors(section_text, section_text1_offset, section_text2_offset, relationship_obj, highlight_dict, section_count);
    console.log(section_text);
    $(this).html(section_text);
    section_count++;
  });

  file_key = relationship_obj.relation_type;
  concepts = {
    'c1': {'text': relationship_obj.c1_text, 'type': relationship_obj.c1_stype},
    'c2': {'text': relationship_obj.c2_text, 'type': relationship_obj.c2_stype}
  };

  var coll = new RelationList(data[file_key]);
  Tree['start'].show( new RelationCompositeView({
    collection: coll,
    concepts: concepts,
    choice: false
  }));

});
