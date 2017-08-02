var resyn_lookup_kinds = {
  'g': 'gene',
  'd': 'disease',
  'c': 'drug'
}

var channel = Backbone.Radio.channel('resyn');


/*
 *  Models & Collections
 */


RESynopsisExtraction = Backbone.Model.extend({
  defaults:{
    'id': null,
    'document_id': null,
    'kind': '',
    'concept_a': {},
    'concept_b': {},
    'answers': [],

    'focus': false,
  }
});


RESynopsisExtractionCollection = Backbone.Collection.extend({
  model: RESynopsisExtraction,
  url: function() {
    return '/task/re/916/analysis/';
  }
});


/*
 * Views
 */


RESynopsisExtractionDetail = Backbone.Marionette.View.extend({
  template: '#re-synopsis-detail-template',

  ui: {
    'context': '#chart-context',
  },

  onRender: function() {
    var self = this;

    var answers = this.model.get('answers');
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
        label = self.model.get('concept_a')['text'] + label + resyn_lookup_kinds[self.model.get('kind'[0])] + ' concept';
      } else if(answer_key == c_2_broken) {
        label = self.model.get('concept_b')['text'] + label + resyn_lookup_kinds[self.model.get('kind')[2]] + ' concept';
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

    console.log(data);

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

  }
});


RESynopsisExtractionItem = Backbone.Marionette.View.extend({
  template: _.template('&#8226;'),

  tagName: 'li',
  className: 'list-inline-item',

  events : {
    'mouseenter': 'triggerDisplay',
    'mousedown': 'triggerDisplay'
  },

  initialize: function(){
    this.listenTo(this.model, 'change:focus', this.render);
  },

  onRender: function() {
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
    if(answers.length) {
      var user_answer_id = _.findWhere(answers, {self: true })['answer']['id'];
      var agree = _.filter(answers, function(obj) { return obj['answer']['id'] == user_answer_id }).length
      var score = agree / answers.length;
      var color = '#45BF55';
      if(score <= .5) { color = '#FFE11A'; }
      if(score <= .25) { color = '#B9121B'; }
      this.$el.css({'color': color});
    }

    if(this.model.get('focus') == true) {
      this.$el.css({'borderBottom': '4px solid gray'});
    } else {
      this.$el.css({'borderBottom': 'none'});
    }

  },

  triggerDisplay: function() {
    channel.trigger('re:synopsis:detail:show', this.model);
  },

});


RESynopsisExtractionView = Backbone.Marionette.CollectionView.extend({
  /* Parent list for RESynopsisExtractionItem */
  tagName: 'ul',
  className: 'list-unstyled list-inline',
  childView: RESynopsisExtractionItem,
});


RESynopsis = Backbone.Marionette.View.extend({
  template: '#re-synopsis-template',
  className: 'row',

  regions: {
    're-extractions': '#re-synopsis-extractions',
    're-extraction-detail': '#re-synopsis-detail-view'
  },

  initialize: function() {
    this.collection = new RESynopsisExtractionCollection({});
    this.collection.fetch();

    this.listenTo(channel, 're:synopsis:detail:show', this.showExtractionDetail)
  },

  onRender: function() {
    this.showChildView('re-extractions', new RESynopsisExtractionView({'collection': this.collection}));
  },

  showExtractionDetail: function(extraction_model) {
    this.collection.invoke('set', {'focus': false});
    extraction_model.set('focus', true);
    this.showChildView('re-extraction-detail', new RESynopsisExtractionDetail({'model': extraction_model}));
  }
});

