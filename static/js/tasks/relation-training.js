Tree.addInitializer(function(options) {
  Tree.addRegions({'start': '#tree-insert'});
  Backbone.Radio.DEBUG = true;
  Tree['convoChannel'] = Backbone.Radio.channel('convo');

  Tree['start'].onShow = function(m, c) {
    if(relation_task_settings['fade']) {
      c.$el.css({'opacity': .25});
      c.$el.append('<div id="tree-cover"></div>');
    }

    setTimeout(function() {
      Tree['convoChannel'].trigger('start', '');
    }, 500);

  }

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
  if (relation_type == 'g_c') {
    $('#c1').addClass('gene');
    $('#c1 .not_correct_stype').text('is not a gene concept?');
    $('#c2').addClass('chemical');
    $('#c2 .not_correct_stype').text('is not a drug concept?');

  } else if (relation_type == 'c_d') {
    $('#c1').addClass('chemical');
    $('#c1 .not_correct_stype').text('is not a drug concept?');
    $('#c2').addClass('disease');
    $('#c2 .not_correct_stype').text('is not a disease concept?');

  } else if (relation_type == 'g_d') {
    $('#c1').addClass('gene');
    $('#c1 .not_correct_stype').text('is not a gene concept?');
    $('#c2').addClass('disease');
    $('#c2 .not_correct_stype').text('is not a disease concept?');
  };

  /* Special for training */
  if( current_relationship.get('concepts')['c1']['text'] == 'Citizen Scientists' ) {
    $('#c1 .not_correct_stype').text('is not a group of people?');
  }

  if( current_relationship.get('concepts')['c2']['text'] == 'Astrology' ) {
    $('#c2 .not_correct_stype').text('is not a field of science?');

  if( current_relationship.get('concepts')['c1']['text'] == 'citizen scientist' ) {
    $('#c1 .not_correct_stype').text('is not a helpful person?');
  }  
  if( current_relationship.get('concepts')['c2']['text'] == 'Biomedical research' ) {
    $('#c2 .not_correct_stype').text('is not a field of study?');
  }

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

    } else {
      /* If no other relations to complete, submit document */
      $('#feedback_modal').on('hidden.bs.modal', function (e) {
        $('#task_relation_results').submit();
      });
    }
  }
});



/* Initalize the page by using the task_data */
var collection;
/* Initalize the Application, but load events later */
Tree.start();

/* Onload request all relation tasks to complete */
collection = new RelationTaskCollection(task_data);
// Sort the collection by C1 on initial data load:
var new_collection = _.sortBy(collection.models, function(c) {
  return c.attributes.concepts.c1.text;
});
collection.models = new_collection;
collection.next();

var current_relationship = collection.findWhere({'current': true});
add_relation_classes(current_relationship);


