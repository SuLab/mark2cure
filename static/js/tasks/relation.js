var data = {
  // "concept-1-broken": "asdfasdfsadf",
  // "concept-2-broken": "asdfasadffsadf",

  "chemical_disease_relation_menu": [{
    "id": "8qota4u8hwtcyp65kz9zm0vjyuxwjt12sko084sn",
    "text": "relates to",
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
      "id": "rhkmksv5jh0vn7p47uk3fwdior6mlgaubwh1l6ow",
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
var section_text1 = "";
var section_text2 = "";
var selected_relation = null;

Tree.addInitializer(function(options) {
  selected_relation = null;
  $('#submit_button').addClass('disabled');

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
    if( selected_relation != null ) {
      $('#submit_button').removeClass('disabled');
    };
    console.log("Selected item:", selected_relation);
  });

  /* When the back toggle is selected */
  Tree['convoChannel'].on('back', function(opts) {
    var collection;

    var parentRel = opts['choice'].get('parentRelation');
    if(parentRel) { collection = parentRel.get('children'); }

    if(opts['collection']) { collection = opts['collection']; }

    /* Backup: Go to the top of the stack */
    collection = coll;
    selected_relation = null;
    $('#submit_button').addClass('disabled');

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
    str = str.replace(new RegExp("^(.{" + String(key_offset) + "})(.{" + String(highlight_dict[key][0]['span']) + "})"),   "$1<strong style='background-color:" + String(color_find(relationship_obj['c'+concept_idx+'_stype'])) + "'>$2</strong>");
    }
    catch(err){
      console.log(err, "LOG THE ERROR");
    };
    };
    return str;
};

function get_section_text(section_text, section_count){
  console.log(section_text, "section text");

  if (section_text1 == "" && section_count == 0) {
    section_text1 = section_text;
  };
  if (section_text2 == "" && section_count == 1) {
    section_text2 = section_text;
  };
  if (section_count == 0){
    section_text = section_text1;
  } else if (section_count == 1) {
    section_text = section_text2;
  };

  return section_text;
};

// Makes object containing text locations
var highlight_dict_LARGE = {};
function get_annotation_spans(relationship_obj) {
  concept_list = [relationship_obj.c1_text, relationship_obj.c2_text];
  type_list = [relationship_obj.c1_stype, relationship_obj.c2_stype];
  location_lists = [relationship_obj.c1_locations, relationship_obj.c2_locations]

  for (i = 0; i < concept_list.length; i++) {
    concept_text = concept_list[i];
    concept_stype = type_list[i];
    location_list = location_lists[i];

    for (locations = 0; locations < location_list.length; locations ++){
        highlight_dict_small = {}
        highlight_list_small = [];

        var location = location_list[locations]
        var location = location.split(":")
        highlight_dict_small['index'] = location[0]
        highlight_dict_small['span'] = location[1]
        highlight_dict_small['text'] = concept_text
        highlight_list_small.push(highlight_dict_small);
        highlight_dict_LARGE[highlight_dict_small['index']] = highlight_list_small;
    };
  };
  return highlight_dict_LARGE;
};

function repeat(str, num_repeats) {
    return new Array(num_repeats + 1).join(str);
};

function add_progress_bar(concept_pairs_remaining, concept_pairs_total) {
  concept_pairs_left = concept_pairs_total - concept_pairs_remaining;
  var circle = "&#8226; ";  // circle symbol
  var remaining_concept_circles = repeat(circle, concept_pairs_left);
  var completed_concept_circles = repeat(circle, concept_pairs_remaining);
  var added_text = $('<span />').css('color', '#7F3CFF');
  added_text.html(remaining_concept_circles);
  $('#testie').html(added_text).css("font-size", "40px").append(completed_concept_circles);
};

function pub_specific_info(relationship_obj) {

  // Third (not guaranteed pubtator).  There are either two or sometimes three pubtators.
  $.ajax({
  url: '/document/pubtator/specific/'+ relationship_obj.pub_list[2] +'.json',
  dataType: 'json',
  async: false,
  success: function(pub3_data) {
      pub3_passage0 = pub3_data.collection.document.passage[0]['annotation'];
      pub3_passage1 = pub3_data.collection.document.passage[1]['annotation'];
      $('.section').each(function(el_idx, el) {
        highlight_dict_LARGE = get_annotation_spans(relationship_obj, pub3_passage0, pub3_passage1);
      });
    },
    error: function() {
      console.log("failed to retrieve third/optional Pubtator API");
      }
  });

  $.ajax({
  url: '/document/pubtator/specific/' + relationship_obj.pub_list[1] + '.json',
  dataType: 'json',
  async: false,
  success: function(pub1_data) {
      pub1_passage0 = pub1_data.collection.document.passage[0]['annotation'];
      pub1_passage1 = pub1_data.collection.document.passage[1]['annotation'];
      $('.section').each(function(el_idx, el) {
        highlight_dict_LARGE = get_annotation_spans(relationship_obj, pub1_passage0, pub1_passage1);
      });
    },
    error: function() {
      console.log("failed to retrieve first required Pubtator API");
      }
  });

  $.ajax({
  url: '/document/pubtator/specific/'+ relationship_obj.pub_list[0] +'.json',
  dataType: 'json',
  async: false,
  success: function(pub2_data) {
    pub2_passage0 = pub2_data.collection.document.passage[0]['annotation'];
    pub2_passage1 = pub2_data.collection.document.passage[1]['annotation'];
    section_text1_offset = pub2_data.collection.document.passage[0]['text'].length
    section_text2_offset = pub2_data.collection.document.passage[1]['text'].length

    // reformat each section and use highlight_dict to change colors in the text at the same time.
    section_count = 0;
    $('.section').each(function(el_idx, el) {

      section_text = $(this).html();
      section_text = get_section_text(section_text, section_count);

      highlight_dict_LARGE = get_annotation_spans(relationship_obj);
      section_text = format_text_colors(section_text, section_text1_offset, section_text2_offset, relationship_obj, highlight_dict_LARGE, section_count);
      $(this).html(section_text).css('line-height', '40px');
      section_count++;
    });
  },
  error: function() {
    console.log("failed to retrieve final required Pubtator API");
    }
  });
  };
// Must go to this API first to retrieve info prior to getting Pub Specific APIs
// Therefore, no point to do this asynchronously with pub_specific_info() function
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
  add_progress_bar(concept_pairs_remaining, concept_pairs_total)
});


var counter = 1;
$('#submit_button').on('click', function(evt) {

  relationship_obj = global_data[counter];  //TODO NEED TO LOG THE CORRECT RELATIONSHIP OBJECT
  previous_relationship_for_ajax = global_data[counter-1];
  counter += 1;
  console.log(relationship_obj.pk);


  $.ajax({
    type: 'POST',
    url: '/relation/test/results/',
    data: $.extend({'csrfmiddlewaretoken': csrf_token },
                   {'relation_type': selected_relation},
                   {'relation': previous_relationship_for_ajax.pk },
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

  add_progress_bar(concept_pairs_remaining, concept_pairs_total)

  section_count = 0;
  $('.section').each(function(el_idx, el) {
    highlight_dict_LARGE = {};
    // TODO make into a function
    section_text = $(this).html();
    section_text = get_section_text(section_text, section_count);

    // Build the higlight dictionary to contain all needed highlights for both
    // concepts in the PAIR
    highlight_dict_LARGE = get_annotation_spans(relationship_obj);
    section_text = format_text_colors(section_text, section_text1_offset, section_text2_offset, relationship_obj, highlight_dict_LARGE, section_count);
    $(this).html(section_text);
    section_count++;
  });

  file_key = relationship_obj.relation_type;
  concepts = {
    'c1': {'text': relationship_obj.c1_text, 'type': relationship_obj.c1_stype},
    'c2': {'text': relationship_obj.c2_text, 'type': relationship_obj.c2_stype}
  };
  $(this).addClass('disabled');
  selected_relation = null;

  var coll = new RelationList(data[file_key]);
  Tree['start'].show( new RelationCompositeView({
    collection: coll,
    concepts: concepts,
    choice: false
  }));

});
