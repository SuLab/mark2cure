if(typeof relation_task_settings !== 'undefined') {

  RelDoc = Backbone.RelationalModel.extend({
    defaults:{
      'focus': false,
    }
  });
  Relations = Backbone.Collection.extend({
    model: RelDoc,
    url: '/task/relation/'+relation_task_settings.document_pk+'/analysis/',
  });

  var color_scale = d3.scale.linear()
      .domain([0, .5, 1])
      .range(['red', 'yellow', 'green']);

  RelationItem = Backbone.Marionette.ItemView.extend({
    template: _.template('&#8226;'),
    tagName: 'li',
    className: 'list-inline-item',

    initialize: function(){
      this.listenTo(this.model, 'change:focus', this.render);
    },

    events : {
      'mouseenter': 'triggerDisplay',
      'mousedown': 'triggerDisplay'
    },

    triggerDisplay: function(evt) {
      var self = this;
      this.model.set('focus', true);
      this.model.collection.each(function(m) {
        if(self.model.id != m.id) {
          m.set({'focus': false});
        }
      });
      Synopsis['convoChannel'].trigger('triggerDisplay', {'model': this.model});
    },

    onRender: function(evt) {
      /*
        # Responses, # Agree, Color, %
        1, 1, Green, 100
        2, 2, Green, 100
        2, 1, Yellow, 50
        3, 3, Green, 100
        3, 2, Green, 66.666,
        3, 1, Yellow, 33
        4, 1, Red, 25
        5, 5, Green, 100
        5, 4, Green, 80
        5, 3, Green, 60
        5, 2, Yellow, 40
        5, 1, Red, 20

        Red = 0 - 25
        Yellow = 26 - 50
        Green = 51 - 100
      */
      var answers = this.model.get('answers');

      var responses = answers.length;
      var user_answer_id = _.findWhere(answers, {self: true })['answer']['id'];
      var agree = _.filter(answers, function(obj) { return obj['answer']['id'] == user_answer_id }).length
      var score = agree / responses;
      var color = '#45BF55';
      if(score <= .5) { color = '#FFE11A'; }
      if(score <= .25) { color = '#B9121B'; }

      this.$el.css({'color': color});

      if(this.model.get('focus') == true) {
        this.$el.css({'borderBottom': '4px solid gray'});
      } else {
        this.$el.css({'borderBottom': 'none'});
      }
    }
  });

  var lookup_kinds = {
    'g': 'gene',
    'd': 'disease',
    'c': 'drug'
  }

  SynopsisCompositeView = Backbone.Marionette.CompositeView.extend({
    template: '#relation-synopsis-template',
    childView: RelationItem,
    childViewContainer: 'ul#relation-synopsis-bar',
    ui: {
      'feedback': '#feedback-next-action-area',
      'feedback_list': '#chart-list',
      'context': '#chart-context',
      'concept_a': '#concept-a',
      'concept_b': '#concept-b',
    },
    initialize: function(evt) {
      var self = this;

      Synopsis['convoChannel'].on('triggerDisplay', function(opts) {
        var model = opts.model;
        self.ui.feedback.html('');
        self.ui.feedback_list.html('');

        self.ui.concept_a.html(model.get('concept_a')['text']);
        self.ui.concept_b.html(model.get('concept_b')['text']);
        self.ui.context.show();

        var answers = model.get('answers');
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
            label = model.get('concept_a')['text'] + label + lookup_kinds[model.get('kind'[0])] + ' concept';
          } else if(answer_key == c_2_broken) {
            label = model.get('concept_b')['text'] + label + lookup_kinds[model.get('kind')[2]] + ' concept';
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

        var chart = d3.select('#feedback-next-action-area').style('width', '100%');
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
      });
    },
  });


  Synopsis = new Backbone.Marionette.Application();
  Synopsis.start();
  var coll = new Relations();
  Synopsis.addInitializer(function(options) {
    coll.fetch();
    Synopsis.addRegions({'start': '#relation-synopsis-insert'});
    Backbone.Radio.DEBUG = true;
    Synopsis['convoChannel'] = Backbone.Radio.channel('convo');
  });
  Synopsis['start'].show(new SynopsisCompositeView({'collection': coll}));

}
