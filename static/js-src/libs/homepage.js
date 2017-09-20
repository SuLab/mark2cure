// $(function() {

// });


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
    return ['col-4 col-sm-3 col-md-2 col-lg-1', 'group-item', 'p-1', group_status].join(' ');
  }
});


HomePageQuestListView = Backbone.Marionette.CollectionView.extend({
  childView: HomePageQuestView,
  childViewEventPrefix: 'questlist',
  className: 'row justify-content-center no-gutters',

  filter: function (child, index, collection) {
    /* Only show the Quests that are enabled */
    // return child.get('enabled');
    return child.get('description') && index < 11;
  }
});


HomePageQuestExplorer = Backbone.Marionette.View.extend({
  template: '#homepage-quest-explorer-template',
  className: 'row',

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

  }
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
