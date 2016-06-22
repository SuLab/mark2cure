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
    .domain([0, 1])
    .range(['red', 'green']);

RelationItem = Backbone.Marionette.ItemView.extend({
  template: _.template('&#8226;'),
  tagName: 'li',

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
    var answers = this.model.get('answers');
    var answer_counts = _.countBy( _.map(answers, function(x) { return x['answer']['id']; }) );
    var user_answer_id = _.findWhere(answers, {self: true })['answer']['id'];

    var answer_count_nums = _.sortBy(answer_counts, function(d) { return -d['value'] });
    var sum = _.reduce(answer_count_nums, function(memo, num){ return memo + num; }, 0);

    var arr = Object.keys( answer_counts ).map(function ( key ) { return answer_counts[key] });
    var majority_total_votes = Math.max.apply( null, arr );
    var majority_answer_id = Object.keys(answer_counts).filter(function(x){ return answer_counts[x] == majority_total_votes; })[0];

    var score = .5;
    /* you match majority and enough data to give green */
    if (user_answer_id == majority_answer_id && majority_total_votes/sum >= 0.51){
      score = 1;

    /* you do not match majority and majority did pretty well */
    } else if ( user_answer_id != majority_answer_id && majority_total_votes/sum >= 0.51 ) {
      score = 0;
    }
    this.$el.css({'color': color_scale(score)});

    if(this.model.get('focus') == true) {
      this.$el.css({'borderBottom': '4px solid gray'});
    } else {
      this.$el.css({'borderBottom': 'none'});
    }
  }
});

SynopsisCompositeView = Backbone.Marionette.CompositeView.extend({
  template: '#relation-synopsis-template',
  childView: RelationItem,
  childViewContainer: 'ul#relation-synopsis-bar',
  ui: {
    'feedback': '#feedback-next-action-area',
    'feedback_list': '#chart-list'
  },
  initialize: function(evt) {
    var self = this;
    Synopsis['convoChannel'].on('triggerDisplay', function(opts) {
      var model = opts.model;
      self.ui.feedback.html('');
      self.ui.feedback_list.html('');

      var answers = model.get('answers');
      var answer_counts = _.countBy( _.map(answers, function(x) { return x['answer']['id']; }) );
      var answer_text = {};
      var personal_ann = '';

      _.each(answers, function(a) {
        answer_text[a.answer.id] = a.answer.text;
        if(a['self']) { personal_ann = a.answer.id; }
      });

      var data = [];
      _.each(_.keys(answer_counts), function(answer_key) {
        data.push({
          'id': answer_key,
          'value': answer_counts[answer_key],
          'label': answer_text[answer_key],
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
