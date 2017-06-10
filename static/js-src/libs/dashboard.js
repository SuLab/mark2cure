


NERDashboard = Backbone.Marionette.View.extend({

  template: '',

  regions: {
    'stats': '',
    'list': '',
  },

  initialize: function() {
    /* > /api/groups/ – list available NER groups
     * > /api/quest/41/ - detail specific NER Group
     */
  }

  onRender: function() {

       var draw_dashboard = function(group, quests) {
        $('#group-'+ group.pk).html('');
        var canvas = d3.select('#group-'+ group.pk);
        var available_quests = _.filter(quests, function(item) { return item.enabled && !item.completed });
        var completion_size = _.map(available_quests, function(item) { return item.completions; });

        var completion_scale = d3.scale.linear()
          .domain([_.min(completion_size), _.max(completion_size)])
          .range(['#00CCFF', '#E64C66']);

        var template = _.template($('#relation-item-template').html());

        var attrs = {
          'class': 'quest col-xs-4 col-sm-3 col-md-3 col-lg-2',
        };
        var styles = {
        };
        var quest = canvas.selectAll('.quest').remove();
        var quest = canvas.selectAll('.quest').data(quests);
        quest.enter().append('div')
          .attr(attrs)
          .style(styles)
          .html(function(d, i) {
            return template({
              'd': d,
              'progress': (d.progress.current/d.progress.required)*100,
              'completions_scale': completion_scale(d.completions),
            });
          });
        quest.transition().attr(attrs);
        quest.exit().remove();
      };

      $.ajax({
        'type': 'GET',
        'url': '/api/groups/',
        'success': function(data) {
          _.each(data, function(v) {
            if(v.enabled) {
              var template = _.template("<div class='row m-t-1'><div class='col-xs-12'><h3><a href='/group/<%- stub %>/'><%- name %></a></h3><a href='/group/<%- stub %>'><p class='text-muted'><%- _.str.prune(description, 100) %></p></a></div><div id='group-<%- pk %>' class='col-xs-12 paragraph-box'><p class='quest-loading text-xs-center'>Loading...</p></div></div>");
              $('#group-selection').append(template(v));

              $.ajax({
                'type': 'GET',
                'url': '/api/quest/'+ v.pk +'/',
                'success': function(data) {
                  draw_dashboard(v, data);

                  $('#group-selection .quest').click(function(evt) {
                    var link = $(this).find('a');
                    if(link.length) { location.href=link.attr('href'); }
                  });

                }
              });

            };

          });
        }
      });



  }

});


REDocument = Backbone.RelationalModel.extend({
  defaults: {}
});

REDocumentList = Backbone.Collection.extend({
  model: Document,
  url: function() { return '/api/relationships/?format=json'; },
});

/*
 * Views
 */
REDocumentView = Backbone.Marionette.View.extend({
  template: '#relation-list-item-template',
  tagName: 'a',
  className: 'relation-item-container',

  onRender : function() {
    this.$el.attr('href', '/task/relation/'+ this.model.get('id') +'/')
  },

});

REDocumentCompositeView = Backbone.Marionette.CollectionView.extend({
  template: _.template('<div id="relation-task-list"></div>'),
  childView  : DocumentView,
  childViewContainer: "#relation-task-list",

  initialize : function(options) {
    this.collection.fetch();
  },

});




REDashboardView = Backbone.Marionette.View.extend({

  template: '',

  regions: {
    'stats': '',
    'list': '',
  },

  onRender: funciotn() {
        if($('#document-relation-board').length) {

        var DocumentRelationBoard = Backbone.Marionette.Application.extend({
          region: '#document-relation-board',

          onStart: function() {
            var docs = new DocumentList();


            var main = this.getRegion();
            main.show( new DocumentCompositeView({'collection': docs}) );
          }
        });
        var doc_rel_board_app = new DocumentRelationBoard();
        doc_rel_board_app.start();

      }



  }

});

LeaderBoardView = Backbone.Marionette.View.extend({

  initialize: function() {
    var settings = new LeaderBoardSettings();
    var users = new UserList();
  },

  onRender: function() {
    var main= this.getRegion();
    main.show( new UserCompositeView({'model': settings, 'collection': users}) );
  }

});

Dashboard = Backbone.Marionette.View.extend({
  template: '#dashboard-template',

  regions: {
    'ner': '#dashboard-ner',
    're': '#dashboard-re',
    'leaderboard': '#dashbaord-leaderboard',
  },

  onRender: function() {
      this.showChildView('ner', new NERDashboardView());
      this.showChildView('re', new REDashboardView());
      this.showChildView('leaderboard', new LeaderBoardView());
  }

});
