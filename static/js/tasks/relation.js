var data = {
  // "concept-1-broken": "asdfasdfsadf",
  // "concept-2-broken": "asdfasadffsadf",

  "chemical_disease_relation_menu": [{
    "id": "8qota4u8hwtcyp65kz9zm0vjyuxwjt12sko084sn",
    "text": "relates to disease",
    "desc": "The chemical has some type of relation to the disease.",
    "example": "Evidence indicates that Amblyomin-X could be a promising candidate for cancer therapy.",
    "children": [{
        "id": "72zuw4bgzz5dniepb3rmls23nsltwocbk274c98m",
        "text": "exacerbates",
        "desc": "The chemical worsens or exacerbates the disease in any manner.",
        "example": "St. John's Wort has been known to worsen symptoms of Parkinson's disease.",
      }, {
        "id": "jilhvc5p2cy0atls8659a1fggjvvkmahwuspy2kr",
        "text": "treats",
        "desc": "The chemical improves or treats the disease in any manner.",
        "example": "Treatment with gingerol may provide a novel therapeutic strategy for prion-mediated neurotoxicity",
      }, {
        "id": "yrrb92b8vtjmcagjj4nx43sbj8wey2moqagk9ea5",
        "text": "increases risk of disease",
        "desc": "Administration of the drug increases chances of getting the disease. Considered an 'adverse side effect' of the drug.",
        "example": "Parkinson's disease risk factors include ageing, exposure to toxins and genetic defects.",
      }, {
        "id": "jyiczzhhupcp7cmebl422ax5dxe1jkwuq647oaw2",
        "text": "may cause",
        "desc": "Administration of the chemical may cause or lead to the disease.",
        "example": "43% of mice given the chemical NVP-BHG712 developed pulmonary metastases.",
      }, {
        "id": "lt18qfxd1ehj7ymxb29wrv6qa41mocwe6eor9dna",
        "text": "prevents",
        "desc": "Administration of the chemical prevents the disease",
        "example": "Malarone is taken daily by travelers wanting to prevent malaria.",
      }]
    }, {
    "id": "4mzrh5ub3nla6b1ostx7qdparjl3lrd9o567ubif",
    "text": "has no relation to",
    "desc": "There is no relation between the chemical and disease.",
    "example": "To study soft tissue sarcoma, an experimental group of mice was given the chemical NVP-BHG712.",
  }],

  "gene_disease_relation_menu": [{
    "id": "qq84lkjfh46gmx4a9n1jpwxwrmbajsy1qctb9u8j",
    "text": "relates to",
    "desc": "The gene is related to the disease in some manner.",
    "example": "The gene responsible for Triple A syndrome, AAAS, has recently been identified.",
    "children": [{
      "id": "u0q779rcrevnu6aki694dqka4fnfwvwqgpl06ybl",
      "text": "altered expression is associated with",
      "desc": "Gene expression is the process by which information from a gene is used in the synthesis of a functional gene product (mRNA, protein, or microRNAs). When gene expression is altered, gene expression is either increased or decreased.",
      "example": "Several studies revealed significantly higher EPHB4 expression in malignancies such as prostate cancer.",
      }, {
      "id": "04110gzdcxz8niuv83ev08ut7lv0xep4iym5sxm5",
      "text": "mutation is associated with",
      "desc": "Gene mutations or aberrations are permanent changes in the gene DNA sequence, differing from the sequence found in most people. Mutations range in size; they can affect anywhere from a single DNA base pair to a large piece of a chromosome that includes multiple genes.",
      "example": "Mutations in the COL5A or COL3A genes are only a few of the genetic causes of Ehlers-Danlos syndrome.",
      }, {
      "id": "5mgmlk7rnbbd8q6amuapd3elwjpnbho0raegv59c",
      "text": "and post-translational modifications are associated with",
      "desc": "Translation is protein synthesis. When changes are made to a protein during or after translation, it is considered a post-translational modification.",
      "example": "A small interfering RNA causes knockdown of ATP2C1 expression, resulting in defects in both post-translational processing of wild-type thyroglobulin",
    }]}, {
    "id": "7aqs9bklotxhbq3r5dcofvsskiefb1yn2nkt1y4a",
    "text": "has no relation to",
    "desc": "The gene does not relate to the disease.",
    "example": "The precise role of Ngly1 in the ERAD process remains unclear in mammals.",
  }],

  "gene_chemical_relation_menu": [{
    "id": "txh8mu2mrnrffik893gr5h0ir7b1y7plgw94n4j7",
    "text": "relates to",

    "children": [{
      "id": "am1wc2yvdcvwcb3yi298xesplbdktzku6wis49iw",
      "text": "amount is changed by",
    }, {
      "id": "5ex6vuro19zeneiwlc8yze6dsq1coxvlpojolwgy",
      "text": "activity is changed by",
    }]

    }, {
    "id": "5mgmlk7rnbbd8q6amuapd3elwjpnbho0raegv59c",
    "text": "has no relation to"
    }
  ]

};

selected_relation = null;

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
    selected_relation = obj['choice'].id;
    console.log("Selected item:", obj['choice'].id);
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

  offset = 0;
  keys = Object.keys(highlight_dict).sort(function(a, b) {return a - b;});

  // process keys in reverse order
  for (i = keys.length - 1; i > -1; i--) {
    key = keys[i];
    key_offset = key;
    if (section_count == 0 && key > section1_length) {  // Passed the end of section 1
      continue;
    }
    else if (section_count == 1 && key < section1_length) {  // In section 2 but key is in section 1
      continue;
    }
    else if (section_count == 1 && key > section1_length) {  // In section 2 and need to shift left by section1_length
      key_offset = key - section1_length - 1;
    }
    else {  // Within section 1
      key_offset = key;
    };
    try {
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
var relationship_obj;
var concepts;
var file_key;
var pub1_passage0, pub1_passage1,
    pub2_passage0, pub2_passage1,
    pub3_passage0, pub3_passage1,
    section_text1_offset, section_text2_offset;


// Makes object containing text locations
var highlight_dict_LARGE = {};
function get_annotation_spans(relationship_obj, passage0, passage1) {
  concept_list = [relationship_obj.c1_text, relationship_obj.c2_text];

  for (i = 0; i < concept_list.length; i++) {
    concept = concept_list[i];

    function find_anns_from_passage(concept, passage){
      highlight_list_small = [];
      for (ann in passage) {
        highlight_dict_small = {};

        if (passage['text'] == concept){
          console.log(passage['text'])
          highlight_dict_small['index'] = passage['location']['@offset'];
          highlight_dict_small['span'] = passage['location']['@length'];
          highlight_dict_small['text'] = passage['text'];
          highlight_list_small.push(highlight_dict_small);
          highlight_dict_LARGE[highlight_dict_small['index']] = highlight_list_small
        } else if (passage[ann]['text'] == concept) {
          console.log(passage[ann]['text'])
          highlight_dict_small['index'] = passage[ann]['location']['@offset'];
          highlight_dict_small['span'] = passage[ann]['location']['@length'];
          highlight_dict_small['text'] = passage[ann]['text'];
          highlight_list_small.push(highlight_dict_small);
          highlight_dict_LARGE[highlight_dict_small['index']] = highlight_list_small;
        };
      };
      return highlight_dict_LARGE;
    };
    highlight_dict_LARGE = find_anns_from_passage(concept, passage0)
    highlight_dict_LARGE = find_anns_from_passage(concept, passage1)


  };
  return highlight_dict_LARGE;
};


function pub_specific_info(relationship_obj) {

    // Third (not guaranteed pubtator).  There are either two or three pubtators.
    $.getJSON('/document/pubtator/specific/'+ relationship_obj.pub_list[2] +'.json', function( pub3_data ) {
      pub3_passage0 = pub3_data.collection.document.passage[0]['annotation'];
      pub3_passage1 = pub3_data.collection.document.passage[1]['annotation'];
      $('.section').each(function(el_idx, el) {
        highlight_dict_LARGE = get_annotation_spans(relationship_obj, pub3_passage0, pub3_passage1);
      });
    });

    // Below are the first and second guaranteed Pubtators. The "third pubtator" is not guaranteed
    $.getJSON('/document/pubtator/specific/'+ relationship_obj.pub_list[1] +'.json', function( pub1_data ) {
      pub1_passage0 = pub1_data.collection.document.passage[0]['annotation'];
      pub1_passage1 = pub1_data.collection.document.passage[1]['annotation'];
      $('.section').each(function(el_idx, el) {
        highlight_dict_LARGE = get_annotation_spans(relationship_obj, pub1_passage0, pub1_passage1);
      });
    });

    $.getJSON('/document/pubtator/specific/'+ relationship_obj.pub_list[0] +'.json', function( pub2_data ) {
      pub2_passage0 = pub2_data.collection.document.passage[0]['annotation'];
      pub2_passage1 = pub2_data.collection.document.passage[1]['annotation'];
      section_text1_offset = pub2_data.collection.document.passage[0]['text'].length
      section_text2_offset = pub2_data.collection.document.passage[1]['text'].length

      // reformat each section and use highlight_dict to change colors in the text at the same time.
      section_count = 0;
      $('.section').each(function(el_idx, el) {
        section_text = $(this).text();
        highlight_dict_LARGE = get_annotation_spans(relationship_obj, pub2_passage0, pub2_passage1);
        section_text = format_text_colors(section_text, section_text1_offset, section_text2_offset, relationship_obj, highlight_dict_LARGE, section_count);
        $(this).html(section_text);
        section_count++;
      });
    });
  };


$.getJSON('/relation/'+ document_pk +'/api/', function(data) {
  global_data = data;

  concept_pairs_remaining = data.length;
  concept_pairs_total = concept_pairs_remaining;

  // API is built off of the "unanswered relations". Therefore this is always the first unanswered one,
  // even if user leaves and returns.
  relationship_obj = data[0];
  pub_specific_info(relationship_obj);

  // Start the tree
  file_key = relationship_obj.relation_type;
  concepts = {
    'c1': {'text': relationship_obj.c1_text, 'type': relationship_obj.c1_stype },
    'c2': {'text': relationship_obj.c2_text, 'type': relationship_obj.c2_stype }
  };
  Tree.start();

});

// Variable praise section:
function yays_selection() {
  yays = ['Sweet!', 'Great!', 'Nice!', 'Awesome!', 'Nice job!', 'Excellent!', 'Wow!', 'Woohoo!', 'Hooray!', 'Yeaahh!', 'Look at you!'];
  yay_word = yays[Math.floor(Math.random()*yays.length)];
  return yay_word;
};

function carry_on_selection(concept_pairs_total, concept_pairs_remaining) {
  carry_on = ['Carry on.', 'Keep going.', "Help some more.", "We still need you.", "Don't stop now.", "", "", "" ];
  if (concept_pairs_total >= 5 && concept_pairs_remaining <= 3 )
    carry_on = ["You're a trooper.", 'Most impressive.', "You're dedicated.", 'Very persistent.', 'Almost there.', 'On a role.', 'You are amazing.'];
  carry_on_word = carry_on[Math.floor(Math.random()*carry_on.length)];
  return carry_on_word;
};

function html_praise_display(concept_pairs_total, concept_pairs_remaining) {
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
  relationship_obj = global_data[counter];
  counter += 1;

  $.ajax({
    type: 'POST',
    url: '/relation/test/results/',
    data: $.extend({'csrfmiddlewaretoken': csrf_token },
                   {'relation_type': "test answer"},
                   {'relation': relationship_obj.pk },
                   {'username': username }),
    cache: false,
    async: false,
    success: function() { console.log('win'); },
    error: function() { console.log('fail'); },
  });

  concept_pairs_remaining -= 1;
  if (concept_pairs_remaining == 0) {
    pass;
  }
  else {
    evt.preventDefault();
  };

  $('#praise_after_completion').hide();
  html_display = html_praise_display(concept_pairs_total, concept_pairs_remaining);
  $('#praise_after_completion').html( html_display ).fadeIn(900);

  section_count = 0;
  $('.section').each(function(el_idx, el) {
    highlight_dict_LARGE = {};
    section_text = $(this).text();

    // Build the higlight dictionary to contain all needed highlights
    highlight_dict_LARGE = get_annotation_spans(relationship_obj, pub1_passage0, pub1_passage1);
    highlight_dict_LARGE = get_annotation_spans(relationship_obj, pub2_passage0, pub2_passage1);
    highlight_dict_LARGE = get_annotation_spans(relationship_obj, pub3_passage0, pub3_passage1);
    console.log(highlight_dict_LARGE, "lg2")
    section_text = format_text_colors(section_text, section_text1_offset, section_text2_offset, relationship_obj, highlight_dict_LARGE, section_count);
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
