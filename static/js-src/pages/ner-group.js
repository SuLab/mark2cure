//
// Models + Collections
//

NERGroupStats = Backbone.Model.extend({
  defaults: {
    "pk": null,
    "enabled": true,
    "stub": "bipg",
    "name": "BiP part two",
    "description": "",

    "total_contributors": 0,
    "document_count": 0,

    "percentage_complete": 0, // What is this?
    "complete_percent": 0,
    "current_avg_f_score": 0,

    "start_date": "",
    'end_date': "in progress",
  },

  url: function() {
    return '/api/ner/list/'+this.group_pk+'/';
  },

  initialize: function(options) {
    this.group_pk = options['group_pk'];
  }
});


NERGroupContributor = Backbone.Model.extend({
  defaults: {
    'username': '',
    'count': 0
  }
});
NERGroupContributorsList = Backbone.Collection.extend({
  model: NERGroupContributor,
  url: '/api/ner/list/31/contributors/'
})

//
// Views
//

NERGroupPageWordCloudView = Backbone.Marionette.View.extend({
  onRender: function() {
  }
});

NERGroupPageNetworkView = Backbone.Marionette.View.extend({
  template: '#ner-group-network',
  className: 'col-12',

  events: {
    'mousedown h4': function() {
      var self = this
      if( $('#network-row').is(":visible")  ) {
        $('#group-network h4 i').removeClass('fa-caret-up').addClass('fa-caret-down');
      } else {
        $('#group-network h4 i').removeClass('fa-caret-down').addClass('fa-caret-up');
      };

      $('#network-row').toggle(function() {
        self.options['sigma'].refresh();
        self.options['sigma'].refresh();
      });
    },

    'mousedown i.fa-plus-circle': function() {
      var c = this.options['sigma'].camera;
      sigma.misc.animation.camera(c, {
        ratio: c.ratio / c.settings('zoomingRatio')
      }, { duration: 250 });
    },

    'mousedown i.fa-minus-circle': function() {
      var c = this.options['sigma'].camera;
      sigma.misc.animation.camera(c, {
        ratio: c.ratio * c.settings('zoomingRatio')
      }, { duration: 250 });
    },

    'mousedown i.fa-rotate-right': function() {
      var c = this.options['sigma'].camera;
      c.goTo({
        angle: c.angle -= .1
      });
    }

  },

  onAttach: function() {
    var self = this;

    this.options['sigma'] = new sigma({
      container: 'network',
      settings: {
        defaultLabelColor: "#000",
        defaultLabelSize: 12,
        defaultLabelBGColor: "#ddd",
        defaultHoverLabelBGColor: "#002147",
        defaultLabelHoverColor: "#fff",

        edgeColor: 'default',
        defaultEdgeColor: '#e7e7e7',

        labelThreshold: 10,
        defaultEdgeType: 'curve',

        hoverFontStyle: "bold",
        fontStyle: "regular",
        activeFontStyle: "regular",

        minNodeSize: 4,
        maxNodeSize: 8,
        minEdgeSize: 1,
        maxEdgeSize: 1,

        zoomMin: .00001,
        zoomMax: 1,
      }
    });
    this.options['filter'] = new sigma.plugins.filter(this.getOption('sigma'));

    $.ajax({
      'type': 'GET',
      'url': '/api/network/'+ self.getOption('group_pk') +'/',
      'success': function(graph) {
        self.options['sigma'].graph.clear();
        self.options['sigma'].graph.read(graph);

        self.options['sigma'].refresh();

        maxDegree = 0;
        self.options['sigma'].graph.nodes().forEach(function(n) {
          maxDegree = Math.max(maxDegree, self.options['sigma'].graph.degree(n.id));
        });
        $('#min-degree').attr('max', maxDegree/3)
      },
      'error': function(d) {
        $('#group-network').html("<div class='row'><div class='col-xs-12'><h4 class='text-xs-center'>Network Unavailable <i class='fa fa-exclamation-triangle fa-1'></i></h4></div></div>");
      }
    });

    this.options['sigma'].refresh();
    this.options['sigma'].bind('hoverNode clickNode', function(e) {
      if(e.type == 'clickNode') {
        var node_color = e.data.node.color;
        var second_search = '';
        if (node_color == '#B1FFA8') {
            second_search = 'gene';
        } else if (node_color == '#d1f3ff') {
            second_search = 'disease';
        } else if (node_color == '#ffd1dc') {
            second_search = 'drug';
        } else {
            second_search = '';
        };
        window.open('https://www.google.com/#safe=off&q='+e.data.node.label+'+'+second_search,'_blank');
      }
    });

    // $('#min-degree').on('change', function() {
    //   var min_degree_num = $(this).val();
    //   $('#min-degree-val').html(min_degree_num);
    //   filter
    //     .undo('min-degree')
    //     .nodesBy(function(n) {
    //       return this.degree(n.id) >= min_degree_num;
    //     }, 'min-degree').apply();
    // });
    // $('#reset-btn').on('click', function() {
    //   filter.undo().apply();
    // });

  }

});

NERGroupPageContributorView = Backbone.Marionette.View.extend({
  template: '#ner-group-homepage-contributor-item-template',
  tagName: 'a',
  className: 'list-group-item justify-content-between',

  onRender : function() {
    this.$el.attr('href', '/u/'+this.model.get('username')+'/');
  }

});

NERGroupPageContributorListView = Backbone.Marionette.CollectionView.extend({
  childView: NERGroupPageContributorView,
  className: 'list-group',
  tagName: 'ul',

  initialize: function() {
    this.collection = new NERGroupContributorsList();
    this.collection.fetch();
  }
});

NERGroupStatsView = Backbone.Marionette.View.extend({
  template: '#ner-group-homepage-stats-template',

  modelEvents: {
    'sync': function() { this.render(); }
  },
});

NERGroupPageView = Backbone.Marionette.View.extend({
  /* this.model NERGroupStats (information about the group)
  *  this.collection NERQuestTaskList (quests within the group) */
  el: '#ner-group-home',
  template: false,

  regions: {
    'stats': '#group-statistics',
    'contributors': '#group-contributors',
    'network': '#group-network-action-area',
    'quest-list': '#group-quest-list'
  },

  initialize: function() {
    this.model = new NERGroupStats(this.options);
    this.model.fetch();

    this.collection = new NERQuestTaskList({'pk': this.getOption('group_pk')});
    this.collection.fetch();
  },

  onRender: function() {
    this.showChildView('stats', new NERGroupStatsView({'model': this.model}));
    this.showChildView('contributors', new NERGroupPageContributorListView());
    this.showChildView('network', new NERGroupPageNetworkView(this.options));
    this.showChildView('quest-list', new NERQuestTaskListView({'collection': this.collection}));
  }
});


