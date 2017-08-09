/*
 * Models + Collections
 */

NERQuestTask = Backbone.Model.extend({
  defaults: {
    'id': 0,
    'name': '1',
    'documents': [],
    'points': 5000,
    'requires_qualification': 7,
    'provides_qualification': 7,
    'meta_url': null,
    'user': {
      'completed': false,
      'enabled': true
    },
    'progress': {
      'current': 15,
      'completed': true,
      'required': 15
    }
  }
});

NERQuestTaskList = Backbone.Collection.extend({
  model: NERQuestTask,
  initialize: function(options) {
    this.pk = options.pk;
    this.fetch();
  },

  url: function() {
    return '/api/ner/list/'+this.pk+'/quests/';
  }
});


NERQuest = Backbone.Model.extend({
  defaults: {},

  // relations: [{
  //   type: 'HasOne',
  //   key: 'details',
  //   relatedModel: NERQuestDetails,
  //   reverseRelation : {
  //     key : 'quest',
  //     type: 'HasOne',
  //     includeInJSON: false,
  //   }
  // }],

  initialize: function() {
    if(this.get('enabled')) {
      this.set('tasks', new NERQuestTaskList({'pk': this.get('pk')}) );
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
  url: '/api/task/stats/',
  defaults: {
    'ner': 0,
    're': 0
  }
});


/*
 * Views
 */

NERQuestTaskView = Backbone.Marionette.View.extend({
  template: '#dashboard-ner-quest-task-template',
  className: 'dashboard-ner-quest-task-item col-3 pb-3',

  templateContext: function() {
    return {
      'percentage_progress': (this.model.get('progress').current/this.model.get('progress').required)*100
    }
  },

});

NERQuestTaskListView = Backbone.Marionette.CollectionView.extend({
  childView: NERQuestTaskView,
  className: 'dashboard-ner-task-list-box row',
  childViewEventPrefix: 'task',
});

NERQuestView = Backbone.Marionette.View.extend({
  template: '#dashboard-ner-quest-template',
  className: 'row dashboard-ner-quest justify-content-center py-4',

  regions: {
    'list': '.task-list-view'
  },

  onRender: function() {
    this.showChildView('list', new NERQuestTaskListView({'collection': this.model.get('tasks')}) );
  }
});


NERQuestListView = Backbone.Marionette.CollectionView.extend({
  childView: NERQuestView,
  className: 'col-12',
  childViewEventPrefix: 'questlist',

  filter: function (child, index, collection) {
    /* Only show the Quests that are enabled */
    return child.get('enabled');
  }
});


NEREmptyListView = Backbone.Marionette.View.extend({
  template: _.template('Sorry, no annotaiotn quests are currently available for you to work on.')
});

NERDashboardUnlockView = Backbone.Marionette.View.extend({
  template: '#dashboard-ner-unlock-template',
  className: 'row justify-content-center',
})

NERDashboardView = Backbone.Marionette.View.extend({
  template: '#dashboard-ner-template',
  className: 'row',

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
  }
});


REDocumentView = Backbone.Marionette.View.extend({
  template: '#relation-list-item-template',
  tagName: 'a',
  className: 'list-group-item list-group-item-action flex-column align-items-start',

  onRender : function() {
    this.$el.attr('href', '/task/re/'+ this.model.get('id') +'/')
  }
});

REDocumentListView = Backbone.Marionette.CollectionView.extend({
  childView: REDocumentView,
  className: 'list-group col-12 my-4',
  childViewEventPrefix: 'relist',

  // filter: function (child, index, collection) {
  //   #<{(| Only show the Quests that are enabled |)}>#
  //   return child.get('enabled');
  // }
});

REEmptyListView = Backbone.Marionette.View.extend({
  template: _.template('Sorry, no relationship tasks currently available for you to work on.')
});

REDashboardUnlockView = Backbone.Marionette.View.extend({
  template: '#dashboard-re-unlock-template',
  className: 'row justify-content-center'
})


REDashboardView = Backbone.Marionette.View.extend({
  template: '#dashboard-re-template',
  className: 'col',

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
    } else {
      this.showChildView('ner', new NERDashboardUnlockView({'model': this.model}) );
    }

    if(this.model.get('re') >= 3) {
      this.showChildView('re', new REDashboardView());
    } else {
      this.showChildView('re', new REDashboardUnlockView({'model': this.model}) );
    }

    this.showChildView('leaderboard', new LeaderBoardView());
  }
});
