/*
 * Models + Collections
 */

NERQuestDetails = Backbone.RelationalModel.extend({
  initialize: function() {
    this.fetch();
  },

  url: function() {
    return '/api/ner/list/'+this.get('quest').get('pk')+'/';
  }
});


NERQuest = Backbone.RelationalModel.extend({
  defaults: {},

  relations: [{
    type: 'HasOne',
    key: 'details',
    relatedModel: NERQuestDetails,
    reverseRelation : {
      key : 'quest',
      type: 'HasOne',
      includeInJSON: false,
    }
  }],

  initialize: function() {
    if(this.get('enabled')) {
      this.set('details', new NERQuestDetails({'quest': this}) );
    }
  }
});


NERQuestList = Backbone.Collection.extend({
  model: NERQuest,
  url: '/api/ner/list/'
});


NERStats = Backbone.Model.extend({
  url: '/api/ner/stats/',
  defaults: {
    'total_score': 0,
    'papers_reviewed': 0,
    'annotations': 0,
    'quests_completed': 0
  }
});


REDocument = Backbone.RelationalModel.extend({
  defaults: {}
});


REDocumentList = Backbone.Collection.extend({
  model: REDocument,
  url: '/api/re/list/',
});


REStats = Backbone.Model.extend({
  url: '/api/re/stats/',
  defaults: {
    'total_score': 0,
    'annotations': 0,
    'quests_completed': 0
  }
});

DashboardTaskStats = Backbone.Model.extend({
  url: '/api/task/stats/'
});


/*
 * Views
 */

NERQuestTaskView = Backbone.Marionette.View.extend({
  template: '#dashboard-ner-quest-task-template',
});

NERQuestTaskListView = Backbone.Marionette.CollectionView.extend({
  childView: NERQuestTaskView,
  childViewEventPrefix: 'task',
});

NERQuestView = Backbone.Marionette.View.extend({
  template: '#dashboard-ner-quest-template',

  regions: {
    'list': '.paragraph-box'
  },

  onRender: function() {
    this.showChildView('list', new NERQuestTaskListView({'collection': this.model.get('details').get('tasks')}) );
  }
});


NERQuestListView = Backbone.Marionette.CollectionView.extend({
  childView: NERQuestView,
  childViewEventPrefix: 'questlist',

  filter: function (child, index, collection) {
    /* Only show the Quests that are enabled */
    return child.get('enabled');
  }
});


NEREmptyListView = Backbone.Marionette.View.extend({
  template: _.template('Sorry, no annotaiotn quests are currently available for you to work on.')
});


NERDashboardView = Backbone.Marionette.View.extend({
  template: '#dashboard-ner-template',
  className: 'col-xs-5',

  regions: {
    'list': '#dashboard-ner-list'
  },

  initialize: function() {
    /* > /api/ner/list/ – list available NER groups
     * > /api/quest/41/ - detail specific NER Group
     */
    this.collection = new NERQuestList();
    this.collection.fetch();
    this.model = new NERStats();
    this.model.fetch();
  },

  modelEvents: {
    'sync': function() { this.render(); }
  },

  onRender: function() {
    this.showChildView('list', new NERQuestListView({'collection': this.collection}));

      //   var available_quests = _.filter(quests, function(item) { return item.enabled && !item.completed });
      //   var completion_size = _.map(available_quests, function(item) { return item.completions; });
      //   var completion_scale = d3.scale.linear()
      //     .domain([_.min(completion_size), _.max(completion_size)])
      //     .range(['#00CCFF', '#E64C66']);
      //   var template = _.template($('#relation-item-template').html());
      //   var attrs = {
      //     'class': 'quest col-xs-4 col-sm-3 col-md-3 col-lg-2',
      //   };
      //   var styles = {
      //   };
      //   var quest = canvas.selectAll('.quest').remove();
      //   var quest = canvas.selectAll('.quest').data(quests);
      //   quest.enter().append('div')
      //     .attr(attrs)
      //     .style(styles)
      //     .html(function(d, i) {
      //       return template({
      //         'd': d,
      //         'progress': (d.progress.current/d.progress.required)*100,
      //         'completions_scale': completion_scale(d.completions),
      //       });
      //     });
      //   quest.transition().attr(attrs);
      //   quest.exit().remove();
      // };
  }
});


REDocumentView = Backbone.Marionette.View.extend({
  template: '#relation-list-item-template',
  tagName: 'a',
  className: 'relation-item-container',

  onRender : function() {
    this.$el.attr('href', '/task/relation/'+ this.model.get('id') +'/')
  }
});

REDocumentListView = Backbone.Marionette.CollectionView.extend({
  childView: REDocumentView,
  childViewEventPrefix: 'relist',

  // filter: function (child, index, collection) {
  //   #<{(| Only show the Quests that are enabled |)}>#
  //   return child.get('enabled');
  // }
});

REEmptyListView = Backbone.Marionette.View.extend({
  template: _.template('Sorry, no relationship tasks currently available for you to work on.')
});


REDashboardView = Backbone.Marionette.View.extend({
  template: '#dashboard-re-template',
  className: 'col-xs-5',

  regions: {
    'list': '#dashboard-re-list',
  },

  modelEvents: {
    'sync': function() { this.render(); }
  },

  initialize: function() {
    /* this.collection =
     * /api/re/list/
     */
    this.collection = new REDocumentList();
    this.collection.fetch();
    this.model = new REStats();
    this.model.fetch();
  },

  onRender: function() {
    this.showChildView('list', new REDocumentListView({'collection': this.collection}));

    // var doc_rel_board_app = new DocumentRelationBoard();
    // doc_rel_board_app.start();
  }
});


DashboardView = Backbone.Marionette.View.extend({
  template: '#dashboard-template',
  className: 'row',

  regions: {
    'ner': '#dashboard-ner',
    're': '#dashboard-re',
    'leaderboard': '#dashbaord-leaderboard',
  },

  initialize: function() {
    this.model = new DashboardTaskStats();
    this.model.fetch();
  },

  modelEvents: {
    'sync': function() { this.render(); }
  },

  onRender: function() {
    if(this.model.get('ner') >= 7) {
      this.showChildView('ner', new NERDashboardView());
    }

    if(this.model.get('re') >= 3) {
      this.showChildView('re', new REDashboardView());
    }

    this.showChildView('leaderboard', new LeaderBoardView());
  }
});
