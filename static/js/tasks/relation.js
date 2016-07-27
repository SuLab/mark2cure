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

var add_relation_classes = function(current_relationship) {
  if (!current_relationship) { return; }

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
          if(annotation) {
            return _.any(annotation.infon, function(infon) {
              return infon['@key'] == 'uid' && _.contains(concept_uids, infon['#text']);
              /* var match = _.filter(concept_uids, function(s) { return infon['#text'].indexOf(s) !== -1 || s.indexOf(infon['#text']) !== -1; }).length;
               * return infon['@key'] == 'uid' && match; */
            });
          } else { return []; }

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
      /* If no other relations to complete, submit document
      * This condition should be accounted for earlier, but
      * keeping here just incase */
      $('#tree-action-area').hide();
      $('#task_relation_submit_document_set').submit();
    }
  }
});


var ProgressItem = Backbone.Marionette.ItemView.extend({
  tagName: 'li',
  className: 'list-inline-item uncompleted',
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


var submit_status = function() {
  /* Submit button status changes accordingly */
  if( Tree.start.currentView.options.choice ) {
    $('#submit_button').attr('disabled', false).removeClass('disabled');
  } else {
    $('#submit_button').attr('disabled', true).addClass('disabled');
  };
};

var lookup_kinds = {
  'g': 'gene',
  'd': 'disease',
  'c': 'drug'
}

function show_results(document_pk, relation_pk) {
  $('#tree-action-area').hide();

  $.getJSON('/task/relation/'+ document_pk +'/analysis/' + relation_pk + '/', function(api_data) {
    var answers = api_data[0]['answers']

    /* obj, key = identifier, value = count */
    var answer_counts = _.countBy( _.map(answers, function(x) { return x['answer']['id']; }) );

    var answer_text = {};
    var personal_ann = '';

    var c_1_broken = "zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy";
    var c_2_broken = "RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer";

    _.each(answers, function(a) {
      answer_text[a.answer.id] = a.answer.text;
      if(a['self']) { personal_ann = a.answer.id; }
    });

    var data = [];
    _.each(_.keys(answer_counts), function(answer_key) {

      var label = answer_text[answer_key];
      if(answer_key == c_1_broken) {
        label = api_data[0]['concept_a']['text'] + label + lookup_kinds[api_data[0]['kind'][0]] + ' concept';
      } else if(answer_key == c_2_broken) {
        label = api_data[0]['concept_b']['text'] + label + lookup_kinds[api_data[0]['kind'][2]] + ' concept';
      }

      data.push({
        'id': answer_key,
        'value': answer_counts[answer_key],
        'label': label,
        'self': answer_key == personal_ann
      });
    });

    data = _.sortBy(data, function(d) { return -d['value'] });

    var max = d3.sum(data.map(function(i){ return i['value']; }));
    var color = d3.scale.category20c();

    var chart = d3.select('#chart').style('width', '100%');
    var bar = chart.selectAll('div')
      .data(data)
      .enter()
        .append('div')
          .attr('class', 'bar-component')
          .style('width', function(d) { return ((d['value']/max)*100) + '%'; } )
          .style('background-color', function(d, i) { return color(i); });

    var bold = function(d, i) {
      if( i == 1 && d['self'] ) {
        return '<strong>';
      } else if( i == 2 && d['self'] ) {
        return '</strong>';
      } else {
        return '';
      }
    }

    var list = d3.select('#chart-list');
    var bar = list.selectAll('li')
      .data(data)
      .enter()
        .append('li')
          .html(function(d, i) { return '<div class="box" style="background-color:'+color(i)+';"></div> <p>' + bold(d, 1) + ((d['value']/max)*100).toFixed() + '% – ' + d['label'] + bold(d, 2) + '</p>'; });

    $('#feedback-next-action-area').addClass('shown').slideDown('fast');
  });
}

$('#next_button').on('click', function(evt) {
  /* Causes the next to load */
  var current_relationship = collection.findWhere({'current': true});
  var next_relationship = collection.findWhere({user_completed: false, available: true});

  if (next_relationship) {
    current_relationship.set('user_completed', true);
    submit_status();
    $('#feedback-next-action-area').hide();
    $('#tree-action-area').show();
    $('#chart').empty();
    $('#chart-list').empty();
  } else {
    $('#tree-action-area').hide();
    $('#task_relation_submit_document_set').submit();
  }

});

var incorrect_flag_arr = ['zl4RlTGwZM9Ud3CCXpU2VZa7eQVnJj0MdbsRBMGy', 'RdKIrcaEOnM4DRk25g5jAfeNC6HSpsFZaiIPqZer'];
$('#submit_button').on('click', function(evt) {
  /* When a user is confirming their Tree choice selection */

  var current_selection = Tree.start.currentView.options.choice;
  var current_relationship = collection.findWhere({'current': true});

  if(current_selection.get('id')) {
    /* Submit the data from the current selection to the server */
    $.ajax({
      type: 'POST',
      url: '/task/relation/'+ relation_task_settings.document_pk +'/'+ current_relationship.id +'/submit/',
      data: $.extend({'csrfmiddlewaretoken': relation_task_settings.csrf_token },
                     {'relation': current_selection.get('id')}),
      cache: false,
      success: function() {

        /* If they flag a concept as incorrect, remove it from the
         * availability pool of tasks for this relationship document */
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

        /* Have them review the community's answers */
        show_results(current_relationship.get('document'), current_relationship.get('id'));
        update_score();

      },
      error: function() { alert('Please refresh your browser and try again.') },
    });
  }
});
