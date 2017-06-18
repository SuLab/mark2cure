$(function() {
  $('a[href*=#]:not([href=#])').click(function() {
    if (location.pathname.replace(/^\//,'') == this.pathname.replace(/^\//,'') && location.hostname == this.hostname) {
      var target = $(this.hash);
      target = target.length ? target : $('[name=' + this.hash.slice(1) +']');
      if (target.length) {
        $('html,body').animate({
          scrollTop: target.offset().top
        }, 1000);
        return false;
      }
    }
  });
});


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


// HomePageStats = Backbone.Model.extend({
//   default: {
//     'annotations': 12345
//   },
//   url: '/api/stats/'
// });


HomePageQuestView = Backbone.Marionette.View.extend({
  template: '#homepage-quest-explorer-item-template',

  onRender: function() {
    console.log('quest item', this.model.attributes);
  }
});


HomePageQuestListView = Backbone.Marionette.CollectionView.extend({
  childView: HomePageQuestView,
  childViewEventPrefix: 'questlist',

  filter: function (child, index, collection) {
    /* Only show the Quests that are enabled */
    return child.get('enabled');
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


// HomePageStatsView = Backbone.Marionette.View.extend({
//   el: '.cta-bar',
//   template: false,
//
//   initialize: function() {
//     this.model = new HomePageStats()
//   },
//
//   onRender: function() {
//     this.showChildView('stats', new HomePageStats());
//     this.showChildView('quests', new HomePageQuestExplorer());
//   }
// });


HomePageView = Backbone.Marionette.View.extend({
  el: '#homepage',
  template: false,

  regions: {
    // 'stats': '#stats-div',
    'quests': '#quest-explorer',
  },

  onRender: function() {
    // this.showChildView('stats', new HomePageStats());
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
