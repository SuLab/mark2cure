Tree.addInitializer(function(options) {
  console.log('Tree INIT:', options);

  Tree.addRegions({'start': '#tree-insert'});

  /* When the app is first loaded */
  var coll = new RelationList(options['data'][options.relation_type]);

  Tree['start'].show( new RelationCompositeView({
    collection: coll,
    concepts: options.concepts,
    choice: false
  }));

  Backbone.Radio.DEBUG = true;
  Tree['convoChannel'] = Backbone.Radio.channel('convo');

  /* When an item is selected */
  Tree['convoChannel'].on('click', function(obj) {
    Tree['start'].show( new RelationCompositeView({
      collection: obj['collection'],
      concepts: options.concepts,
      choice: obj['choice']
    }));
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
      concepts: options.concepts,
      choice: false
    }));
  });
});


RelationTask = Backbone.Model.extend({
  defaults: {
      id: null,
      document: null,
      relation_type: null,
      concepts: {},
      user_completed: false
    }
});

RelationTaskCollection = Backbone.Collection.extend({
  model: RelationTask,
  initialize: function (models, options) {
    this.on('change:user_completed', this.next, this);
  },

  next: function() {
    var next_relationship = collection.findWhere({user_completed: false});
    if (next_relationship) {

      this.current_relationship = collection.findWhere({user_completed: false});
      var coll = new RelationList( relation_task_settings['data'][ this.current_relationship.get('relation_type') ] );
      Tree['start'].show( new RelationCompositeView({
        collection: coll,
        concepts: this.current_relationship.get('concepts'),
        choice: false
      }));

    } else {
      $('#task_relation_results').submit();
    }
  }

  //remaining
});


var collection;
$.getJSON('/task/relation/'+ relation_task_settings.document_pk +'/api/', function(data) {
  if(relation_task_settings) {
    var filtered_data = _.where( data, {'user_completed': false} );
    collection = new RelationTaskCollection(filtered_data);

    var current_relationship = collection.findWhere({user_completed: false});
    collection.current_relationship = current_relationship;
    Tree.start( $.extend(collection.current_relationship.toJSON(), relation_task_settings) );

  } else {
    alert('Sorry, there was a problem.');
  }
});


$('#submit_button').on('click', function(evt) {
  var current_selection = Tree.start.currentView.options.choice.get('id');

  if(current_selection) {
    $.ajax({
      type: 'POST',
      url: '/task/relation/'+ relation_task_settings.document_pk +'/'+ collection.current_relationship.id +'/submit/',
      data: $.extend({'csrfmiddlewaretoken': relation_task_settings.csrf_token },
                     {'relation': current_selection}),
      cache: false,
      success: function() {
        collection.current_relationship.set('user_completed', true);
      },
      error: function() { alert('Please refresh your browswer and try again.') },
    });
  }


});
