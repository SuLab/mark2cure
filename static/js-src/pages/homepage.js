var sci_comm_data = {
  "nodes": [
    {
      "id": "n0",
      "label": "Sudan",
      "size": 3,
      "x": -3,
      "y": -3,
      "color": "#0cf"
    },
    {
      "id": "n1",
      "label": "Uganda",
      "size": 3,
      "x": -3,
      "y": -2,
      "color": "#0cf"
    },
    {
      "id": "n2",
      "label": "Democratic Republic of the Congo",
      "size": 3,
      "x": -2,
      "y": -2,
      "color": "#0cf"
    },
    {
      "id": "n3",
      "label": "Central African Republic",
      "size": 3,
      "x": -3,
      "y": 3,
      "color": "#7F3CFF"
    },
    {
      "id": "n4",
      "label": "Congo",
      "size": 3,
      "x": -2,
      "y": -2,
      "color": "#7F3CFF"
    },
    {
      "id": "n5",
      "label": "Gabon",
      "size": 3,
      "x": -1,
      "y": -1,
      "color": "#7F3CFF"
    },
    {
      "id": "n6",
      "label": "Ape Diseases",
      "size": 3,
      "x": 0,
      "y": 0,
      "color": "#333"
    },
    {
      "id": "n7",
      "label": "Cote d'Ivoire",
      "size": 3,
      "x": 0,
      "y": -1,
      "color": "#E64C66"
    },
    {
      "id": "n8",
      "label": "Monkey Diseases",
      "size": 3,
      "x": 1,
      "y": 1,
      "color": "#333"
    },
    {
      "id": "n9",
      "label": "Vaccines, DNA",
      "size": 3,
      "x": 2,
      "y": -4,
      "color": "#333"
    },
    {
      "id": "n10",
      "label": "Immunity, Humoral",
      "size": 3,
      "x": 3,
      "y": -2,
      "color": "#333"
    },
  ],
  "edges": [
    {
      "id": "e0",
      "source": "n0",
      "target": "n2"
    },
    {
      "id": "e1",
      "source": "n1",
      "target": "n2"
    },
    {
      "id": "e2",
      "source": "n2",
      "target": "n5"
    },
    {
      "id": "e3",
      "source": "n3",
      "target": "n4"
    },
    {
      "id": "e4",
      "source": "n4",
      "target": "n5"
    },
    {
      "id": "e5",
      "source": "n5",
      "target": "n6"
    },
    {
      "id": "e6",
      "source": "n6",
      "target": "n7"
    },
    {
      "id": "e7",
      "source": "n6",
      "target": "n7"
    },
    {
      "id": "e8",
      "source": "n6",
      "target": "n8"
    },
    {
      "id": "e9",
      "source": "n8",
      "target": "n9"
    },
    {
      "id": "e10",
      "source": "n9",
      "target": "n10"
    }
  ]
};

HomePageQuest = Backbone.Model.extend({
  defaults: {
    'pk': null,
    'name': '',
    'stub': '',
    'description': '',
    'enabled': false,
    'complete_percent': 0.0
  },
});


HomePageQuestList = Backbone.Collection.extend({
  model: HomePageQuest,
  url: '/api/ner/list/'
});


Mark2CureStats = Backbone.Model.extend({
  default: {
    'ner_annotations': 0,
    're_annotations': 0
  },
  url: '/api/mark2cure/stats/',

  numberWithCommas: function(x) {
    return Math.round(x).toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  }
});


HomePageQuestView = Backbone.Marionette.View.extend({
  template: '#homepage-quest-explorer-item-template',
  className: function() {
    var group_status = this.model.get('enabled') ? 'active-group' : 'ended-group';
    return ['col-10 col-sm-5 col-md-3 col-lg-2 col-xl-1', 'group-item', 'p-1', group_status].join(' ');
  }
});


HomePageQuestListView = Backbone.Marionette.CollectionView.extend({
  childView: HomePageQuestView,
  childViewEventPrefix: 'questlist',
  className: 'row justify-content-center no-gutters',

  filter: function (child, index, collection) {
    /* Only show the Quests that are enabled */
    // return child.get('enabled');
    return child.get('description') && index < 10;
  }
});


HomePageQuestExplorer = Backbone.Marionette.View.extend({
  template: '#homepage-quest-explorer-template',
  className: 'row justify-content-center',

  regions: {
    'list': '.list'
  },

  initialize: function() {
    this.collection = new HomePageQuestList();
    this.collection.fetch();
  },

  onRender: function() {
    this.showChildView('list', new HomePageQuestListView({'collection': this.collection}));
  }
});


HomePageNetwork = Backbone.Marionette.View.extend({
  template: '#homepage-network',
  className: 'col-12',

  onAttach: function() {
    var self = this;

    var s = new sigma({
      graph: sci_comm_data,
      container: 'network',
      settings: {
        // defaultNodeColor: '#ec5148',
        edgeColor: 'default',
        defaultEdgeColor: '#e7e7e7',

        labelSizeRatio: .5,

        touchEnabled: false,
        mouseEnabled: false,
        mouseWheelEnabled: false,
        doubleClickEnabled: false,
        eventsEnabled: false,
        enableHovering: false

        // defaultLabelColor: "#000",
        // defaultLabelSize: 12,
        // defaultLabelBGColor: "#ddd",
        // defaultHoverLabelBGColor: "#002147",
        // defaultLabelHoverColor: "#fff",
        // labelThreshold: 10,
        // defaultEdgeType: 'curve',
        // hoverFontStyle: "bold",
        // fontStyle: "regular",
        // activeFontStyle: "regular",
        // minNodeSize: 4,
        // maxNodeSize: 8,
        // minEdgeSize: 1,
        // maxEdgeSize: 1,
      }
    });

    s.startForceAtlas2({
      worker: true,
      barnesHutOptimize: false,
      scalingRatio: 2,
      gravity: .125,
      outboundAttractionDistribution: true,
    });

    setTimeout(function() {
      s.killForceAtlas2();
    }, 5000);
  }

});

HomePageView = Backbone.Marionette.View.extend({
  el: '#homepage',
  template: false,

  ui: {
    'annotations': '#annotation-counter',
    'watch_more': '#watch-more-link'
  },

  events: {
    'mousedown @ui.watch_more': function(evt) {
      var target = $('#watch-more');
      $('html,body').animate({
        scrollTop: target.offset().top
      }, 1000);
    }
  },

  regions: {
    'quests': '#quest-explorer',
    'network': '#landing-network-browser'
  },

  initialize: function() {
    this.model = new Mark2CureStats();
    this.model.fetch();
  },

  modelEvents: {
    'sync': function() {
      var total_annotations = this.model.get('ner_annotations') + this.model.get('re_annotations');
      this.ui.annotations.text( this.model.numberWithCommas(total_annotations) );
    }
  },

  onRender: function() {
    this.showChildView('quests', new HomePageQuestExplorer());
    this.showChildView('network', new HomePageNetwork());
  },



});


/* ETC (find a better place for this) */

$(function() {
  var footerHeight = $('.footer').height();
  $('.out').css('margin-bottom', -footerHeight);
  $('.push').css('height', footerHeight);
});

$(window).resize(function() {
  var footerHeight = $('.footer').height();
  $('.out').css('margin-bottom', -footerHeight);
  $('.push').css('height', footerHeight);
});
