Tree.addInitializer(function(options) {
  Tree.addRegions({'start': '#tree-insert'});
  Backbone.Radio.DEBUG = true;
  Tree['convoChannel'] = Backbone.Radio.channel('convo');
});


RelationTask = Backbone.Model.extend({
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

// submit button status changes accordingly
var submit_status = function() {
  if( Tree.start.currentView.options.choice ) {
    $('#submit_button').attr('disabled', false).removeClass('disabled');
  } else {
    $('#submit_button').attr('disabled', true).addClass('disabled');
  };
};

var add_relation_classes = function(current_relationship) {
  var relation_type = current_relationship.get('relation_type');
  if (relation_type == "g_c") {
    $('#c1').addClass('gene');
    $('#c1 .not_correct_stype').text('is not a gene concept?');
    $('#c2').addClass('chemical');
    $('#c2 .not_correct_stype').text('is not a drug concept?');

  } else if (relation_type == "c_d") {
    $('#c1').addClass('chemical');
    $('#c1 .not_correct_stype').text('is not a drug concept?');
    $('#c2').addClass('disease');
    $('#c2 .not_correct_stype').text('is not a disease concept?');

  } else if (relation_type == "g_d") {
    $('#c1').addClass('gene');
    $('#c1 .not_correct_stype').text('is not a gene concept?');
    $('#c2').addClass('disease');
    $('#c2 .not_correct_stype').text('is not a disease concept?');
  };
};

RelationTaskCollection = Backbone.Collection.extend({
  model: RelationTask,
  initialize: function (models, options) {
    this.on('change:user_completed', this.next, this);
  },
  next: function() {
    var next_relationship = collection.findWhere({user_completed: false, available: true});
    if (next_relationship) {
      var self = this;
      if( Tree['convoChannel'] ) { Tree['convoChannel'].unbind(); }

      /* Assign an uncompleted relationship as the current focused task */
      collection.each(function(r) { r.set('current', false); })
      next_relationship.set('current', true);
      var current_relationship = collection.findWhere({'current': true});
      var concepts = current_relationship.get('concepts');

      var current_idx = collection.indexOf(current_relationship);
      if(current_idx >= 1) {
        var previous_concepts = collection.at(current_idx-1).get('concepts');

        if(previous_concepts['c1'].id != concepts['c1'].id) { concepts['c1']['fadeIn'] = true; };
        if(previous_concepts['c2'].id != concepts['c2'].id) { concepts['c2']['fadeIn'] = true; };
      }

      var coll = new RelationList( relation_task_settings['data'][ current_relationship.get('relation_type') ] );
      var view_options = {collection: coll, concepts: current_relationship.get('concepts'), choice: false, first_draw: true };

      Tree['start'].show(new RelationCompositeView(view_options));
      add_relation_classes(current_relationship);

      /* When an item is selected */
      Tree['convoChannel'].on('click', function(obj) {
        var rcv = new RelationCompositeView({collection: obj['collection'], concepts: concepts, choice: obj['choice'], first_draw: false });
        Tree['start'].show(rcv);
        submit_status();
        add_relation_classes(current_relationship);
      });

      /* When the back toggle is selected */
      Tree['convoChannel'].on('back', function(opts) {
        /* Clicking back always completely resets. Backup: Go to the top of the stack */
        view_options['first_draw'] = false;
        Tree['start'].show(new RelationCompositeView(view_options));
        submit_status();
        add_relation_classes(current_relationship);
      });

      /* When C1 or C2 is incorrect */
      Tree['convoChannel'].on('error', function(obj) {
        var rcv = new RelationCompositeView({
          collection: new RelationList([]),
          concepts: concepts,
          choice: new Backbone.Model( relation_task_settings['data'][ obj+'_broken' ] ),
          first_draw: false
        })
        Tree['start'].show(rcv);
        submit_status();
        add_relation_classes(current_relationship);
      });


      var concept_uids = [concepts['c1'].id, concepts['c2'].id];
      tmp_passages = [];

      _.each(passages, function(p, p_idx) {
        /* Deep clone passage objects */
        tmp_passage = $.extend({}, p);

        tmp_passage['annotation'] = _.filter(tmp_passage.annotation, function(annotation) {
          return _.any(annotation.infon, function(infon) {
            return infon['@key'] == "uid" && _.contains(concept_uids, infon['#text']);
          });
        });

        var p = new Paragraph({'text': tmp_passage.text});
        YPet[''+p_idx].show( new WordCollectionView({
          collection: p.get('words'),
          passage_json: tmp_passage,
        }) );
        YPet[''+p_idx].currentView.drawBioC(tmp_passage, false);
        YPet[''+p_idx].currentView.drawBioC(null, true);

      });

    } else {
      /* If no other relations to complete, submit document */
      $('#task_relation_results').submit();
    }
  }
});


var ProgressItem = Backbone.Marionette.ItemView.extend({
  tagName: 'li',
  className: 'uncompleted',
  template: _.template('&#8226;'),
  initialize : function(options) {
    this.listenTo(this.model, 'change:user_completed', this.render);
    this.listenTo(this.model, 'change:current', this.render);
    this.listenTo(this.model, 'change:available', this.render);
  },
  onRender: function() {
    var class_options = ['uncompleted', 'completed', 'skip', 'active']
    var self = this;
    _.each(class_options, function(c) { self.$el.removeClass(c); });

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


var ProgressView = Backbone.Marionette.CollectionView.extend({
  childView: ProgressItem
});


/* Initalize the page by loading all relation tasks
 * and fetching all required data */
var collection;
$.getJSON('/task/relation/'+ relation_task_settings.document_pk +'/api/', function(data) {
  if(relation_task_settings) {

    /* Initalize the Application, but load events later */
    Tree.start();

    /* Onload request all relation tasks to complete */
    collection = new RelationTaskCollection(data);
    // Sort the collection by C1 on initial data load:
    var new_collection = _.sortBy(collection.models, function(c) {
      return c.attributes.concepts.c1.text;
    });
    collection.models = new_collection;
    collection.next();

    var current_relationship = collection.findWhere({'current': true});
    add_relation_classes(current_relationship);
    /* Init Progressbar event listner */
    (new ProgressView({
      collection: collection,
      el: '#progress-bar'
    })).render();

  } else {
    /* There were no relation tasks to complete for this document (scoped to this user) */
    alert('Sorry, there was a problem.');
  }
});


var incorrect_flag_arr = ["zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy", "RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer"];
$('#submit_button').on('click', function(evt) {
  var current_selection = Tree.start.currentView.options.choice;
  var current_relationship = collection.findWhere({'current': true});

  if(current_selection.get('id')) {
    $.ajax({
      type: 'POST',
      url: '/task/relation/'+ relation_task_settings.document_pk +'/'+ current_relationship.id +'/submit/',
      data: $.extend({'csrfmiddlewaretoken': relation_task_settings.csrf_token },
                     {'relation': current_selection.get('id')}),
      cache: false,
      success: function() {

        if( _.contains(incorrect_flag_arr, current_selection.get('id')) ) {
          var incorrect_key = 'c'+(1+incorrect_flag_arr.indexOf(current_selection.get('id')));
          var incorrect_concept = current_relationship.get('concepts')[ incorrect_key ];

          collection.each(function(relation) {
            var concepts = relation.get('concepts');
            if( concepts['c1'].id == incorrect_concept['id'] || concepts['c2'].id == incorrect_concept['id'] ) {
              relation.set('available', false);
            }
          });

        }

        current_relationship.set('user_completed', true);
        submit_status();

      },
      error: function() { alert('Please refresh your browswer and try again.') },
    });
  }
});
