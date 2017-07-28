var incorrect_id_arr = ['zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy', 'RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer'];

var channel = Backbone.Radio.channel('tree');

/*
 *  Models & Collections
 */

REChoice = Backbone.RelationalModel.extend({
  /* The comparison between A => B */
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
  /* The relationship tree of options from a
   * pre-defined selection pool */
  model: REChoice,
  url: '/api/v1/words',
});

REConcept = Backbone.Model.extend({
  /* The Concept term that is being verified and compared to another */
  defaults: {
    id: null,
    text: null,
    index: null,
    type: null
  },
  initialize: function(obj) {
    if(obj.type.length<=1) {

      try {
        var val = {
          'g': 'gene',
          'd': 'disease',
          'c': 'drug'
        }[obj.type]
      } catch(err) {
        throw 'Value not an acceptable option'
      }

      this.set('type', val);
    }
  }
});

REExtractionAnswer = Backbone.Model.extend({
  /*
  * The A => B submission answer
  */
  defaults: {
    percentage: 0.0,
    self: false,
    color: '#000'
  }
});

REExtractionAnswerList = Backbone.Collection.extend({
  /* The list of REExtractionAnswer results for a single Relationship ID */
  model: REExtractionAnswer,
  url: '',
  // data = _.sortBy(data, function(d) { return -d['value'] });
  comparator: 'value',

  initialize: function() {
    var self = this;
    this.color = d3.scale.category20c();

    this.listenTo(this, 'add', function(answer) {
      this.max_value = d3.sum( self.pluck('value') );

      this.each(function(a) {
        a.set('percentage', Math.round(((a.get('value')/self.max_value)*100)*10)/10);
        a.set('color', self.color( self.indexOf(a)) )
      });

    });
  },
});

REDocumentResult = Backbone.RelationalModel.extend({
  /* Provides the information for RE Document results */
  defaults: {
    'task': {
      'pk': null,
      'name': ''
    },
    'award': {
      'pk': null,
      'amount': 0
    },
    'uqr_created': false,
    'uqr_pk': null
  }
});

REExtraction = Backbone.Model.extend({
  /*
  * The A => B proposition
  */
  defaults: {
    id: null,
    relation_type: null,
    concepts: {},
    user_completed: false,
    available: true,
    current: false,
    results: null,
  }
});

REExtractionList = Backbone.Collection.extend({
  /* The list of REExtraction events for a single Document */
  model: REExtraction,

  get_active: function() {
    /* Return the next RE Task for a Document */
    var next_relationship = this.findWhere({user_completed: false, available: true});
    if (next_relationship) {
      /* Assign an uncompleted relationship as the current focused task */
      this.each(function(r) { r.set('current', false); })
      next_relationship.set('current', true);
      return next_relationship;
    } else {
      console.log('trigger that!', this);
      this.trigger('re:document:complete');
      return false;
    }
  }
});

/*
 * Views
 */

REProgressItem = Backbone.Marionette.View.extend({
  /* Drop down list of REChoices */
  template: _.template('&#8226;'),
  tagName: 'li',
  className: 'list-inline-item',

  onRender: function() {
    if ( !this.model.get('available') ) {
      this.$el.addClass('skip');
    }
    if ( this.model.get('user_completed') ) {
      this.$el.addClass('completed');
    }
    if ( this.model.get('current') ) {
      this.$el.addClass('active');
    }
  }
});

REProgressView = Backbone.Marionette.CollectionView.extend({
  /* Parent list for REProgressItem */
  tagName: 'ul',
  className: 'list-unstyled list-inline',
  childView: REProgressItem,
  childViewEventPrefix: 'progress'
});

RENavigationView = Backbone.Marionette.View.extend({
  /* The gray progression bar and score indicator
  * - Modal: None
  * - Collection: REExtractionList
  */
  template: '#tree-navigation-template',
  className: 'row justify-content-between',

  regions: {
    'progress': '#progress-bar'
  },

  ui: {
    'score': '#score'
  },

  initialize: function() {
    this.listenTo(channel, 'tree:navigation:score:update', this.updateScore);
  },

  onRender: function() {
    this.showChildView('progress', new REProgressView({'collection': this.collection}));
    this.updateScore();
  },

  updateScore: function() {
    var self = this;

    $.ajax({
      type: 'GET',
      url: '/u/points/',
      headers: {'X-CSRFTOKEN': this.options.csrf_token},
      cache: false,
      async: false,
      success: function(data) {
        self.ui.score.html(data.points);
      }
    });
  }
});

REExtractionAnswerItem = Backbone.Marionette.View.extend({
  template: '#reextraction-result-answer-item-template',
  tagName: 'li',
});

REExtractionAnswerView = Backbone.Marionette.CollectionView.extend({
  tagName: 'ul',
  className: 'list-unstyled',
  childView: REExtractionAnswerItem,
});


RECompletedView = Backbone.Marionette.View.extend({
  /* Landing page for when RE Document is complete (do next, tweet, fb share, etc)
   * this.model = REQuestResult
   * this.collection = None
   */
  template: '#tree-completed-template',
  className: 'row',

  ui: {
    'next': '#next'
  },

  events: {
    'mousedown @ui.next': 'nextQuest'
  },

  initialize: function() {
    this.model = new REDocumentResult({});
    var self = this;


    var self = this;
    $.ajax({
      type: 'POST',
      url: '/task/relation/'+ self.model.get('document_id') +'/',
      data: {'csrfmiddlewaretoken': self.options.csrf_token },
      cache: false,
      success: function(api_data) {
        // self.triggerMethod('reconfirm:confirm:success', {'child_view': childView, 'data': api_data});

        //     self.model = new NERQuestResult(data)
        //     channel.trigger('ypet:navigation:score:update');
        //     self.render();

      },
    });



  },

  nextQuest: function() {
    /* Use API to determine next Quest for User */
    // $.ajax({
    //   type: 'GET',
    //   url: '/api/ner/list/'+ this.model.get('group').pk +'/quests/',
    //   headers: {'X-CSRFTOKEN': this.options.csrf_token},
    //   success: function(data) {
    //     var set = _.filter(data, function(q) {
    //       return q.progress.completed == false && q.user.completed == false && q.user.enabled == true;
    //     });
    //     if(set.length) {
    //       window.location = '/task/entity-recognition/quest/' + set[0].id;
    //     } else { window.location = '/dashboard/'; }
    //   }
    // });

  }
});



REExtractionResultView = Backbone.Marionette.View.extend({
  /* The results interface for a single REExtraction
  * – Modal: REExtraction
  * - Collection: REExtractionAnswerList */
  template: '#reextraction-results-template',
  className: 'row my-3 justify-content-center',
  ui: {
    'answers_chart': 't',
    'button': 'button'
  },

  regions: {
    'answers-list': '#reextraction-answers-list'
  },

  triggers: {
    'mousedown @ui.button': 're:extractionresult:next'
  },

  onAttach: function() {
    var self = this;
    this.showChildView('answers-list', new REExtractionAnswerView({'collection': this.collection}));

    var chart = d3.select('#reextraction-answers-chart').style('width', '100%');
    var bar = chart.selectAll('div')
      .data(this.collection.toJSON())
      .enter()
        .append('div')
          .attr('class', 'bar-component')
          .style('width', function(d) { return ((d['value']/self.collection.max_value)*100) + '%'; } )
          .style('background-color', function(d, i) { return self.collection.color(i); });

  }
});

REConfirmView = Backbone.Marionette.View.extend({
  /* The actionable interface for submitting
  * the selected REChoice.
  *
  * Note: Used as a mechanism for storing the final REChoice to submit
  *
  * – Modal: REChoice
  * - Collection: None */
  template: _.template('<button id="submit_button" class="btn btn-primary btn-block disabled" disabled="disabled">Submit</button>'),
  tagName: 'div',
  className: 'col-10 col-md-4',
  ui: {
    'button': 'button'
  },
  triggers: {
    'mousedown @ui.button': 'reconfirm:confirm'
  },
  onRender: function() {
    if(this.model) {
      this.ui.button.attr('disabled', false).removeClass('disabled');
    }
  }
});

RESelectedChoiceView = Backbone.Marionette.View.extend({
  /* The display for showing the currently
  * selected REChoice, or a prompt to take action
  * – Modal: REChoice
  * - Collection: None */
  template: _.template(''),
  tagName: 'h3',
  triggers : {
    'mousedown': 'reselectedchoice:clear'
  },
  onRender: function() {
    if(this.model) {
      this.$el.html(this.model.get('text'));
      this.$el.attr('class', 'relation-go-back');
    } else {
      this.$el.html('Select a Relationship below...');
      this.$el.attr('class', 'disabled');
    }
  }
});

REChoiceView = Backbone.Marionette.View.extend({
  /* Drop down list of REChoices */
  template: _.template('<%= text %>'),
  tagName: 'a',
  className: 'list-group-item',

  triggers : {
    'mousedown': 'click'
  }
});

REChoicesView = Backbone.Marionette.CollectionView.extend({
  /* Parent list for REChoices */
  tagName: 'ul',
  childView: REChoiceView,
  childViewEventPrefix: 'rechoice'
});

RELoadingView = Backbone.Marionette.View.extend({
  /* Initial HTML before a REExtractionList is available
  * - Model = None
  * - Collection = None
  */
  template: _.template('<div class="loader"></div>'),
  className: 'loader-container'
});

REConceptView = Backbone.Marionette.View.extend({
  /* The display + actionable area for the
  * REConcepts being compared
  * - Modal: REConcept
  * - Collection: None
  */
  template: '#reextraction-concept-template',
  className: 'concept row',

  ui: {
    'flag': 'div.flag'
  },
  events: {
    'mouseover @ui.flag': function(evt) {
      this.$el.addClass('incorrect');
    },
    'mouseout @ui.flag': function(evt) {
      this.$el.removeClass('incorrect');
    }
  },
  triggers: {
    'mousedown .flag': 'reconcept:incorrect',
  },

  onRender: function() {
    this.$el.addClass(this.model.get('type'));

    /* Determine if the concept names shoudl fade in (to indicate that they're new) */
    // var current_relationship = collection.findWhere({'current': true});
    // var concepts = current_relationship.get('concepts');
    // var current_idx = collection.indexOf(current_relationship);
    // if(current_idx >= 1) {
    //   var previous_concepts = collection.at(current_idx-1).get('concepts');
    //   if(previous_concepts['c1'].id != concepts['c1'].id) { concepts['c1']['fadeIn'] = true; };
    //   if(previous_concepts['c2'].id != concepts['c2'].id) { concepts['c2']['fadeIn'] = true; };
    // }
    // if(this.options.first_draw) {
    //   if( _.has(concepts['c1'], 'fadeIn') ) {
    //     this.ui.c1.addClass('fade-in one');
    //   };
    //   if( _.has(concepts['c2'], 'fadeIn') ) {
    //     this.ui.c2.addClass('fade-in one');
    //   };
    // }
  }
});

REExtractionView = Backbone.Marionette.View.extend({
  /* The tool for describing how
  * A relates to B
  * this.model = a REExtraction instance
  * this.collection = REChoices
  */
  template: '#reextraction-template',
  className: 'row justify-content-center',

  regions: {
    'c1': '#c1',
    'selected_choice': '#selected-choice',
    'list': '#rechoices-list',
    'c2': '#c2',
    'confirm': '#tree-confirm'
  },

  childViewEvents: {
    'rechoice:click': function(childView, evt) {
      /* When a REChoice is selected
       * 1. Set the RESelectionView (backbutton) text
       * 2. Set the list updated to any REChoice children
       * 3. Set the confirmation to the REChoice
       */
      this.showChildView('selected_choice', new RESelectedChoiceView({'model': childView.model}));
      this.showChildView('list', new REChoicesView({'collection': childView.model.get('children')}));
      this.showChildView('confirm', new REConfirmView({'model': childView.model}));
    },
    'reselectedchoice:clear': function(childView, evt) {
      /* Clicking back always completely resets (top of the stack) the RESelection */
      this.showChildView('selected_choice', new RESelectedChoiceView({'model': null}));
      var choices_collection = new REChoices(relation_data[this.model.get('relation_type')]);
      this.showChildView('list', new REChoicesView({'collection': choices_collection}));
      this.showChildView('confirm', new REConfirmView({'model': null}));
    },
    'reconcept:incorrect': function(childView, evt) {
      /* Flag a concept as being incorrect
      * 1. Display informative text in the RESelectionView
      * 2. Clear any potential REChoice listings
      * 3. Set the incorrect concept REChoice to the confirmation
      */
      var m = new REChoice(relation_data['c_'+childView.model.get('index')+1+'_broken'])
      m.set('text', childView.model.get('text') + ' is not a ' + childView.model.get('type') + ' concept');
      this.showChildView('selected_choice', new RESelectedChoiceView({'model': m}));
      this.showChildView('list', new REChoicesView({'collection': new REChoices([])}));
      this.showChildView('confirm', new REConfirmView({'model': new REChoice(relation_data['c_'+(childView.model.get('index')+1)+'_broken'])}));
    },
    'reconfirm:confirm': function(childView, evt) {
      /* Submit the selected REChoice to the server */
      var self = this;
      $.ajax({
        type: 'POST',
        url: '/task/relation/'+ self.model.get('document_id') +'/'+ self.model.get('id') +'/submit/',
        data: $.extend({'csrfmiddlewaretoken': self.options.csrf_token }, {'relation': childView.model.get('id')}),
        cache: false,
        success: function(api_data) {
          self.triggerMethod('reconfirm:confirm:success', {'child_view': childView, 'data': api_data});
          channel.trigger('tree:navigation:score:update');
        },
        error: function() { self.triggerMethod('reconfirm:confirm:error'); },
      });
    },
  },

  onRender: function() {
    /* Prepare the REExtraction subviews
    * 1. Initalize Concept 1 and Concept 2
    * 2. Initalize an empty RESelectionView
    * 3. Initalize the REChoices for this comparison type
    * 4. Initalize the Confirmation with an empty value (how we store the selected answer)
    */
    var c_a = new REConcept(_.extend(this.model.get('concepts')['c1'], {index: 0}));
    this.showChildView('c1', new REConceptView({'model': c_a}));

    var c_b = new REConcept(_.extend(this.model.get('concepts')['c2'], {index: 1}));
    this.showChildView('c2', new REConceptView({'model': c_b}));

    this.showChildView('selected_choice', new RESelectedChoiceView({'model': null}));
    var choices_collection = new REChoices(relation_data[this.model.get('relation_type')]);
    this.showChildView('list', new REChoicesView({'collection': choices_collection}));
    this.showChildView('confirm', new REConfirmView({'model': null}));
  }
});

Tree = Backbone.Marionette.View.extend({
  /* The top level view for all interations of
  * - Model = current / active REExtraction event
  * - Collection = REExtractionList
  */
  template: '#tree-template',
  className: 'row justify-content-center',

  regions: {
    'navigation': '#tree-navigation',
    'extraction': '#tree-selection',
    'extraction-results': '#tree-selection-results',
    'text': '#tree-text',
    // 'footer': '#footer'
  },

  childViewEvents: {
    'reconfirm:confirm:success': function(obj) {
      /* Selected REChoice submitted to the server
      * X 1. Start to download the results
      * 2. Flag future REExtraction models of concept corrections
      * 3. Show REExtraction results
      */
      var self = this;
      var answers_collection = new REExtractionAnswerList([]);

      // Extend any naming labels
      var answer_text = {};
      _.each(obj.data[0]['answers'], function(a) {
        var val = a.answer.text;
        if( _.contains(incorrect_id_arr, a.answer.id) ) {
          var flagged_concept = new REConcept( self.model.get('concepts')['c'+(_.indexOf(incorrect_id_arr, a.answer.id)+1)] );
          val = flagged_concept.get('text')+' is not a '+ flagged_concept.get('type') +' concept';
        }
        answer_text[a.answer.id] = val;
      });

      // Add them all to the REExtractionAnswersList
      var answer_counts = _.countBy( _.map(obj.data[0]['answers'], function(x) { return x['answer']['id']; }) );
      _.each(_.keys(answer_counts), function(answer_key) {
         answers_collection.add(new REExtractionAnswer({
           'id': answer_key,
           'value': answer_counts[answer_key],
           'label': answer_text[answer_key],
           'self': answer_key == obj.child_view.model.get('id')
        }) );
      });
      answers_collection.sort();

      /* Submitting results for a REExtraction */
      this.getRegion('extraction').empty();
      this.showChildView('navigation', new RENavigationView({'collection': this.collection}));
      this.showChildView('extraction-results', new REExtractionResultView({'model': this.model, 'collection': answers_collection}));
    },
    'reconfirm:confirm:error': function(childView) {
      alert('error error');
    },
    're:extractionresult:next': function(childView) {
      /* Go next after reviewing the results */
      this.emptyRegions();

      this.model = this.collection.get_active();
      if(this.model) {
        this.showChildView('navigation', new RENavigationView({'collection': this.collection}));
        this.options['model'] = this.model;
        this.showChildView('extraction', new REExtractionView(this.options));
        this.showChildView('text', new YPet(this.options));
      }
    },
  },

  initialize: function() {
    /* Perform any library state preparation like loading collections
    * or differentiating between Task and Training
    */
    this.options.training = false;
    this.options.mode = 're';
    this.collection = new REExtractionList();

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
        self.model = self.collection.get_active();

        self.render();
      });

    } else {
      this.model = this.collection.get_active();
      this.render();
    }
  },

  collectionEvents: {
    're:document:complete': function() {
      console.log('re : document : complete');
      // When the user just completed the last NER Task
      this.getRegion('extraction').empty();
      this.showChildView('results', new NERQuestCompletedView(this.options));
    },
  },

  onRender: function() {
    if(this.collection && this.model) {
      this.options['model'] = this.model;
      this.options['mode'] = 're';
      this.showChildView('navigation', new RENavigationView({'collection': this.collection}));
      this.showChildView('extraction', new REExtractionView(this.options));
      this.showChildView('text', new YPet(this.options));
    } else {
      this.showChildView('extraction', new RELoadingView());
    };
  }
});
