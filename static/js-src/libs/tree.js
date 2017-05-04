/*
 *  Models & Collections
 */

REChoice = Backbone.RelationalModel.extend({
  defaults: {
    id: '',
    text: '',
    selected: false
  },

  relations: [{
    type: 'HasMany',
    key: 'children',

    relatedModel: 'REChoice',
    collectionType: 'REChoices',

    reverseRelation : {
      key : 'parentREChoice',
      includeInJSON: false,
    }
  }],

});

REChoices = Backbone.Collection.extend({
  model: REChoice,
  url: '/api/v1/words',
});


REExtraction = Backbone.Model.extend({
  /*
  * The A => B proposition
  */
  defaults: {
      id: null,
      document: null,
      relation_type: null,
      concepts: {},
      user_completed: false,
      available: true,
      current: false
    }
});

REExtractionList = Backbone.Collection.extend({
  model: REExtraction,
  // initialize: function (models, options) {
  //   this.on('change:user_completed', this.next, this);
  // },
  next: function() {
    var next_relationship = this.findWhere({user_completed: false, available: true});
    if (next_relationship) {
      /* Assign an uncompleted relationship as the current focused task */
      this.each(function(r) { r.set('current', false); })
      next_relationship.set('current', true);
      return next_relationship;
    } else {
      return false;
    }
  }
});


    // #<{(| Sort the collection by C1 on initial data load |)}>#
    // var new_collection = _.sortBy(collection.models, function(c) {
    //   return c.attributes.concepts.c1.text;
    // });
    // collection.models = new_collection;
    // collection.next();


/*
 * ETC
 */
// var add_relation_classes = function(current_relationship) {
//   if (!current_relationship) { return; }
//
//   var relation_type = current_relationship.get('relation_type');
//   if (relation_type == 'g_c') {
//     $('#c1').addClass('gene');
//     $('#c1 .not_correct_stype').text('is not a gene concept?');
//     $('#c2').addClass('chemical');
//     $('#c2 .not_correct_stype').text('is not a drug concept?');
//
//   } else if (relation_type == 'c_d') {
//     $('#c1').addClass('chemical');
//     $('#c1 .not_correct_stype').text('is not a drug concept?');
//     $('#c2').addClass('disease');
//     $('#c2 .not_correct_stype').text('is not a disease concept?');
//
//   } else if (relation_type == 'g_d') {
//     $('#c1').addClass('gene');
//     $('#c1 .not_correct_stype').text('is not a gene concept?');
//     $('#c2').addClass('disease');
//     $('#c2 .not_correct_stype').text('is not a disease concept?');
//   };
//
//   #<{(| Special for training |)}>#
//   if( current_relationship.get('concepts')['c1']['text'] == 'Citizen Scientists' ) {
//     $('#c1 .not_correct_stype').text('is not a group of people?');
//   }
//   if( current_relationship.get('concepts')['c2']['text'] == 'Astrology' ) {
//     $('#c2 .not_correct_stype').text('is not a field of science?');
//   }
//   if( current_relationship.get('concepts')['c1']['text'] == 'citizen scientist' ) {
//     $('#c1 .not_correct_stype').text('is not a helpful person?');
//   }
//   if( current_relationship.get('concepts')['c2']['text'] == 'Biomedical research' ) {
//     $('#c2 .not_correct_stype').text('is not a field of study?');
//   }
// };
//
// function get_stype(the_current_relation) {
//   return $(this).hasClass('.gene');
// };

// function show_results(document_pk, relation_pk) {
//   $('#tree-action-area').hide();
//
//   $.getJSON('/task/relation/'+ document_pk +'/analysis/' + relation_pk + '/', function(api_data) {
//     var answers = api_data[0]['answers']
//
//     #<{(| obj, key = identifier, value = count |)}>#
//     var answer_counts = _.countBy( _.map(answers, function(x) { return x['answer']['id']; }) );
//
//     var answer_text = {};
//     var personal_ann = '';
//
//     var c_1_broken = "zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy";
//     var c_2_broken = "RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer";
//
//     _.each(answers, function(a) {
//       answer_text[a.answer.id] = a.answer.text;
//       if(a['self']) { personal_ann = a.answer.id; }
//     });
//
//     var data = [];
//     _.each(_.keys(answer_counts), function(answer_key) {
//
//       var label = answer_text[answer_key];
//       if(answer_key == c_1_broken) {
//         label = api_data[0]['concept_a']['text'] + label + lookup_kinds[api_data[0]['kind'][0]] + ' concept';
//       } else if(answer_key == c_2_broken) {
//         label = api_data[0]['concept_b']['text'] + label + lookup_kinds[api_data[0]['kind'][2]] + ' concept';
//       }
//
//       data.push({
//         'id': answer_key,
//         'value': answer_counts[answer_key],
//         'label': label,
//         'self': answer_key == personal_ann
//       });
//     });
//
//     data = _.sortBy(data, function(d) { return -d['value'] });
//
//     var max = d3.sum(data.map(function(i){ return i['value']; }));
//     var color = d3.scale.category20c();
//
//     var chart = d3.select('#chart').style('width', '100%');
//     var bar = chart.selectAll('div')
//       .data(data)
//       .enter()
//         .append('div')
//           .attr('class', 'bar-component')
//           .style('width', function(d) { return ((d['value']/max)*100) + '%'; } )
//           .style('background-color', function(d, i) { return color(i); });
//
//     var bold = function(d, i) {
//       if( i == 1 && d['self'] ) {
//         return '<strong>';
//       } else if( i == 2 && d['self'] ) {
//         return '</strong>';
//       } else {
//         return '';
//       }
//     }
//
//     var list = d3.select('#chart-list');
//     var bar = list.selectAll('li')
//       .data(data)
//       .enter()
//         .append('li')
//           .html(function(d, i) { return '<div class="box" style="background-color:'+color(i)+';"></div> <p>' + bold(d, 1) + ((d['value']/max)*100).toFixed() + '% – ' + d['label'] + bold(d, 2) + '</p>'; });
//
//     $('#feedback-next-action-area').addClass('shown').slideDown('fast');
//   });
// }

var incorrect_id_arr = ['zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy', 'RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer'];
// var lookup_kinds = {
//   'g': 'gene',
//   'd': 'disease',
//   'c': 'drug'
// }

/*
 * Views
 */

// var ProgressItem = Backbone.Marionette.View.extend({
//   tagName: 'li',
//   className: 'list-inline-item uncompleted',
//   template: _.template('&#8226;'),
//   initialize : function(options) {
//     this.listenTo(this.model, 'change:user_completed', this.render);
//     this.listenTo(this.model, 'change:current', this.render);
//     this.listenTo(this.model, 'change:available', this.render);
//   },
//   onRender: function() {
//     var class_options = ['uncompleted', 'completed', 'skip', 'active']
//     var self = this;
//     _.each(class_options, function(c) { self.$el.removeClass(c); });
//
//     if ( !this.model.get('available') ) {
//       this.$el.addClass('skip');
//     }
//     if ( this.model.get('user_completed') ) {
//       this.$el.addClass('completed');
//     }
//     if ( this.model.get('current') ) {
//       this.$el.addClass('active');
//     }
//   }
// });

// var ProgressView = Backbone.Marionette.CollectionView.extend({
//   childView: ProgressItem
// });

REConfirm = Backbone.Marionette.View.extend({
  template: _.template('<button id="submit_button" class="btn btn-primary btn-block disabled" disabled="disabled">Submit</button>'),
  tagName: 'div',
  className: 'col-xs-10 col-xs-offset-1 col-sm-10 col-sm-offset-1 col-md-4 col-md-offset-4',
  ui: {
    'button': 'button'
  },
  triggers: {
    'mousedown @ui.button': 'rechoice:selected:confirm'
  },
  onRender: function() {
    if(this.model) {
      this.ui.button.attr('disabled', false).removeClass('disabled');
    }
  }
});

RESelectedChoiceView = Backbone.Marionette.View.extend({
  template: _.template(''),
  tagName: 'h3',
  triggers : {
    'mousedown': 'rechoice:selected:click'
  },
  onRender: function() {
    if(this.model) {
      this.$el.html(this.model.get('text'));
    } else {
      this.$el.html('Select a Relationship below...');
      this.$el.attr('class', 'disabled');
    }
  }
});

REChoiceView = Backbone.Marionette.View.extend({
  template: _.template('<%= text %>'),
  tagName: 'a',
  className: 'list-group-item',

  triggers : {
    'mousedown': 'click'
  }

});

REChoicesView = Backbone.Marionette.CollectionView.extend({
  tagName: 'ul',
  childView: REChoiceView,
  childViewEventPrefix: 'rechoice'
});


RELoadingView = Backbone.Marionette.View.extend({
  template: _.template('<h1>Loading!</h1>'),
});

REExtractionView = Backbone.Marionette.View.extend({
  /* The tool for describing how
  * A relates to B
  * this.model = a REExtraction instance
  * this.collection = REChoices
  */
  template: '#reextraction-template',

  regions: {
    'list': '#rechoices-list',
    'selected_choice': '#selected-choice',
    'confirm': '#tree-confirm'
  },

  ui: {
    'c1': '#c1',
    'c2': '#c2',
    'c1_not_correct': '#c1_not_correct',
    'c2_not_correct': '#c2_not_correct',
    'incorrect_buttons': '.fa.fa-times-circle',
  },

  events : {
    // 'mouseover @ui.incorrect_buttons': 'hoverIncorrect',
    // 'mouseout @ui.incorrect_buttons': 'hoverOutIncorrect',
    'mousedown @ui.c1_not_correct': 'c1Error',
    'mousedown @ui.c2_not_correct': 'c2Error',
  },

  // hoverIncorrect: function(evt) {
  //   $(evt.target).parent().addClass('incorrect');
  // },
  //
  // hoverOutIncorrect: function(evt) {
  //   var choice = this.options['choice'];
  //   if( !_.contains(incorrect_id_arr, choice.id) ) {
  //     $(evt.target).parent().removeClass('incorrect');
  //   }
  // },

  c1Error: function(evt) {
      /* When C1 is incorrect */
      var model = new REChoice(relation_data.c_1_broken);
      this.showChildView('confirm', new REConfirm({'model': model}));
      this.showChildView('list', new REChoicesView({'collection': new REChoices([])}));
      this.showChildView('selected_choice', new RESelectedChoiceView({'model': model}));
  },

  c2Error: function(evt) {
      /* When C2 is incorrect */
      var model = new REChoice(relation_data.c_2_broken);
      this.showChildView('confirm', new REConfirm({'model': model}));
      this.showChildView('list', new REChoicesView({'collection': new REChoices([])}));
      this.showChildView('selected_choice', new RESelectedChoiceView({'model': model}));
  },

  childViewEvents: {
    'rechoice:click': function(childView, evt) {
      /* When an item is selected */
      /*
       * 1. Set the current choice to the ID
       * 2. Set the backbutton reference if available
       * 3. Set the list updated to any children
       *
       */
      var children = childView.model.get('children');
      this.showChildView('confirm', new REConfirm({'model': childView.model}));
      this.showChildView('selected_choice', new RESelectedChoiceView({'model': childView.model}));
      this.showChildView('list', new REChoicesView({'collection': children}));
    },

    'rechoice:selected:click': function(childView, evt) {
      /* When the back toggle is selected */
      /* Clicking back always completely resets. Backup: Go to the top of the stack */
      var choices_collection = new REChoices(relation_data['g_d']);
      this.showChildView('confirm', new REConfirm({'model': null}));
      this.showChildView('selected_choice', new RESelectedChoiceView({'model': null}));
      this.showChildView('list', new REChoicesView({'collection': choices_collection}));
    },

    'rechoice:selected:confirm': function(childView, evt) {
      var self = this;
      /* When the selection was confirmed, submit the data from the current selection to the server */
      $.ajax({
        type: 'POST',
        url: '/task/relation/'+ self.model.get('document_id') +'/'+ self.model.get('id') +'/submit/',
        data: $.extend({'csrfmiddlewaretoken': self.options.csrf_token },
                         {'relation': childView.model.get('id')}),
        cache: false,
        success: function() {

          /* If they flag a concept as incorrect, remove it from the
           * availability pool of tasks for this relationship document */
          if( _.contains(incorrect_id_arr, childView.model.get('id')) ) {
            self.triggerMethod('reextraction:concept:flagged', childView.model);
          }

          /* Have them review the community's answers */
          // show_results(current_relationship.get('document_id'), current_relationship.get('id'));
          // update_score();
        },
        error: function() { alert('Please refresh your browser and try again.') },
      });


      // var choices_collection = new REChoices(relation_data['g_d']);
      // this.showChildView('confirm', new REConfirm({'model': null}));
      // this.showChildView('selected_choice', new RESelectedChoiceView({'model': null}));
      // this.showChildView('list', new REChoicesView({'collection': choices_collection}));
    }

  },

  onRender: function() {
    var choices_collection = new REChoices(relation_data['g_d']);

    this.showChildView('confirm', new REConfirm({'model': null}));
    this.showChildView('selected_choice', new RESelectedChoiceView({'model': null}));
    this.showChildView('list', new REChoicesView({'collection': choices_collection}));

    /* Determine if the concept names shoudl fade in (to indicate that they're new) */
    // var current_relationship = collection.findWhere({'current': true});
    // var concepts = current_relationship.get('concepts');
    // var current_idx = collection.indexOf(current_relationship);
    // if(current_idx >= 1) {
    //   var previous_concepts = collection.at(current_idx-1).get('concepts');
    //   if(previous_concepts['c1'].id != concepts['c1'].id) { concepts['c1']['fadeIn'] = true; };
    //   if(previous_concepts['c2'].id != concepts['c2'].id) { concepts['c2']['fadeIn'] = true; };
    // }

    // var self = this;
    // var choice = this.options['choice'];
    // var concepts = this.options['concepts'];
    //
    // if(this.options.first_draw) {
    //   if( _.has(concepts['c1'], 'fadeIn') ) {
    //     this.ui.c1.addClass('fade-in one');
    //   };
    //   if( _.has(concepts['c2'], 'fadeIn') ) {
    //     this.ui.c2.addClass('fade-in one');
    //   };
    // }
    //
    // if(choice) {
    //   function get_stype_word(stype){
    //     if (stype == 'g'){
    //       stype_word = 'gene';
    //     } else if ( stype == 'c'){
    //       stype_word = 'drug';
    //     } else {
    //       stype_word = 'disease';
    //     }
    //     return stype_word;
    //   };
    //
    //   if (choice.id == "zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy") {
    //     this.ui.relation.removeClass('disabled').text( concepts['c1'].text + choice.get('text') + get_stype_word(concepts['c1'].type) + ' concept' );
    //     this.ui.c1.addClass('incorrect');
    //
    //   } else if (choice.id == "RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer") {
    //     this.ui.relation.removeClass('disabled').text( concepts['c2'].text + choice.get('text') + get_stype_word(concepts['c2'].type) + ' concept' );
    //     this.ui.c2.addClass('incorrect');
    //
    //     #<{(| Training specific
    //      * (TODO) get out of here
    //      |)}>#
    //     if(concepts['c2']['text'] == 'Astrology') {
    //       this.ui.relation.removeClass('disabled').text( concepts['c2'].text + choice.get('text') + 'field of science');
    //     };
    //
    //   } else {
    //     this.ui.relation.removeClass('disabled').text( choice.get('text') );
    //   };
    // };
    //
    // if(choice || this.collection.parentREChoice) {
    //   this.ui.relation.removeClass('disabled');
    //   this.ui.relation.addClass('relation-go-back');
    // }
  }

});

RelationTextView = Backbone.Marionette.View.extend({
  template: _.template('<p>Text / YPet will go here</p>'),
  onRender: function() {
    // RelationApp['start'].show(new RelationCompositeView(view_options));
    // add_relation_classes(current_relationship);
    //
    // // This to the else is different from the training (it's not in the training)
    // var concept_uids = [concepts['c1'].id, concepts['c2'].id];
    // tmp_passages = [];
    //
    // _.each(passages, function(p, p_idx) {
    //   #<{(| Deep clone passage objects |)}>#
    //   tmp_passage = $.extend({}, p);
    //
    //   tmp_passage['annotation'] = _.filter(tmp_passage.annotation, function(annotation) {
    //     if(annotation) {
    //       return _.any(annotation.infon, function(infon) {
    //         return infon['@key'] == 'uid' && _.contains(concept_uids, infon['#text']);
    //         #<{(| var match = _.filter(concept_uids, function(s) { return infon['#text'].indexOf(s) !== -1 || s.indexOf(infon['#text']) !== -1; }).length;
    //          * return infon['@key'] == 'uid' && match; |)}>#
    //       });
    //     } else { return []; }
    //
    //   });
    //
    //   var p = new Paragraph({'text': tmp_passage.text});
    //   YPet[''+p_idx].show( new WordCollectionView({
    //     collection: p.get('words'),
    //     passage_json: tmp_passage,
    //   }) );
    //   YPet[''+p_idx].currentView.drawBioC(tmp_passage, false);
    //   YPet[''+p_idx].currentView.drawBioC(null, true);
    // });
  }
});

Tree = Backbone.Marionette.View.extend({
  /* The top level view for all interations of
  * the relationship interface
  */
  template: '#tree-template',

  regions: {
    'progress': '#tree-progress',
    'selection': '#tree-selection',
    'selection-results': '#tree-selection-results',
    'text': '#tree-text'
  },

  childViewEvents: {
    'selection:submit': function(childView) {
      /* Submitting results for a REExtraction */
      this.emptyRegions();
      this.showChildView('selection-results', new View({}));
    },

    'selection:results:next': function(childView) {
      /* Go next after reviewing the results */
      this.emptyRegions();
      this.showChildView('selection-results', new REExtractionView({}));

      /* If no other relations to complete, submit document
      * This condition should be accounted for earlier, but
      * keeping here just incase */
      // $('#tree-action-area').hide();
      // $('#task_relation_submit_document_set').submit();

      // From the relation-training.js file
      /* If no other relations to complete, submit document */
      // $('#feedback_modal').on('hidden.bs.modal', function (e) {
      //   $('#task_relation_results').submit();
      // });

    },

    'reextraction:concept:flagged': function(childView) {
      console.log(childView);
      var incorrect_key = 'c'+(1+incorrect_id_arr.indexOf(childView.id));
      console.log('flagged', this, incorrect_key);
      // var incorrect_concept = current_relationship.get('concepts')[ incorrect_key ];

      // collection.each(function(relation) {
      //   var concepts = relation.get('concepts');
      //   if( concepts['c1'].id == incorrect_concept['id'] || concepts['c2'].id == incorrect_concept['id'] ) {
      //     relation.set('available', false);
      //   }
      // });
    }
  },

  initialize: function() {
    this.options.training = false;
    if( _.keys(_.pick(this.options, 'csrf_token', 'document_pk', 'document_pmid')).length != 3) {
      this.options.training = true
    };

    if(!this.options.training && !this.collection) {
      var self = this;
      /* Initalize the page by loading all relation tasks
       * and fetching all required data */
      $.getJSON('/task/relation/'+ this.options.document_pk +'/api/', function(data) {
        /* Onload request all relation tasks to complete */
        self.collection = new REExtractionList(data);

        /* Sort the collection by C1 on initial data load */
        var new_collection = _.sortBy(self.collection.models, function(c) {
          return c.attributes.concepts.c1.text;
        });
        self.collection.models = new_collection;
        self.model = self.collection.next();

        self.render();
      });

    } else {
      this.model = this.collection.next();
    }

  },

  onRender: function() {
    // this.showChildView('progress', new View({'collection': this.collection}));

    if(this.collection && this.model) {
      this.options['model'] = this.model;
      this.showChildView('selection', new REExtractionView(this.options));
      this.showChildView('text', new RelationTextView(this.options));

    } else {
      this.showChildView('selection', new RELoadingView());
    };

  }

});
