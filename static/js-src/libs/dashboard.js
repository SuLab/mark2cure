/*
 * Models + Collections
 */

NERQuest = Backbone.Collection.extend({
});

NERQuestList = Backbone.Collection.extend({
  model: NERQuest,
});

REDocument = Backbone.RelationalModel.extend({
  defaults: {}
});

REDocumentList = Backbone.Collection.extend({
  model: REDocument,
  url: function() { return '/api/relationships/?format=json'; },
});

/*
 * Views
 */
NERQuestView = Backbone.Marionette.View.extend({
  template: _.template('GOD QUEST.')
});


NERQuestList = Backbone.Marionette.CollectionView.extend({
  tagName: 'ul',
  className: '',
  childView: NERQuestView,
  childViewEventPrefix: 'questlist'
});


NEREmptyListView = Backbone.Marionette.View.extend({
  template: _.template('Sorry, no annotaiotn quests are currently available for you to work on.')
});


NERDashboardView = Backbone.Marionette.View.extend({
  template: '#dashboard-ner-template',

  regions: {
    'list': '#dashboard-ner-list'
  },

  initialize: function() {
    /* > /api/groups/ – list available NER groups
     * > /api/quest/41/ - detail specific NER Group
     */
  },

  onRender: function() {
    this.showChildView('list', new NERQuestList());
      //  var draw_dashboard = function(group, quests) {
      //   $('#group-'+ group.pk).html('');
      //   var canvas = d3.select('#group-'+ group.pk);
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
      // $.ajax({
      //   'type': 'GET',
      //   'url': '/api/groups/',
      //   'success': function(data) {
      //     _.each(data, function(v) {
      //       if(v.enabled) {
      //         var template = _.template("<div class='row m-t-1'><div class='col-xs-12'><h3><a href='/group/<%- stub %>/'><%- name %></a></h3><a href='/group/<%- stub %>'><p class='text-muted'><%- _.str.prune(description, 100) %></p></a></div><div id='group-<%- pk %>' class='col-xs-12 paragraph-box'><p class='quest-loading text-xs-center'>Loading...</p></div></div>");
      //         $('#group-selection').append(template(v));
      //         $.ajax({
      //           'type': 'GET',
      //           'url': '/api/quest/'+ v.pk +'/',
      //           'success': function(data) {
      //             draw_dashboard(v, data);
      //             $('#group-selection .quest').click(function(evt) {
      //               var link = $(this).find('a');
      //               if(link.length) { location.href=link.attr('href'); }
      //             });
      //           }
      //         });
      //       };
      //     });
      //   }
      // });
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


REEmptyListView = Backbone.Marionette.View.extend({
  template: _.template('Sorry, no relationship tasks currently available for you to work on.')
});


REDashboardView = Backbone.Marionette.View.extend({
  template: '#dashboard-re-template',

  regions: {
    'stats': '',
    'list': '#dashboard-re-list',
  },

  onRender: function() {
    // var docs = new DocumentList();
    // this.showChildView('list', new DocumentCompositeView({'collection': docs}) );
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

  onRender: function() {
    this.showChildView('ner', new NERDashboardView());
    // this.showChildView('re', new REDashboardView());
    // this.showChildView('leaderboard', new LeaderBoardView());
  }
});
