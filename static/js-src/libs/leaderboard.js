/*
 *  Models & Collections
 */
LeaderBoardSettings = Backbone.RelationalModel.extend({
  defaults: {'time_range': 0, 'text': 'This Month', 'days': 31, 'api': 'users'},
  time_key: [
      {'text': 'This Month', 'days': 31},
      {'text': 'This Week', 'days': 7},
      {'text': 'Today', 'days': 1}
  ],

  toggle_api: function() {
    if(this.get('api') == 'users') {
      this.set('api', 'teams');
    } else {
      this.set('api', 'users');
    }
  },

  toggle_time_range: function() {
    var idx = this.get('time_range');
    idx++;
    if(idx==3) { idx = 0; }
    this.set('time_range', idx);
  },
});

LeaderBoardUser = Backbone.RelationalModel.extend({
  defaults: {'user': false, 'hover': false}
});

LeaderBoardUserList = Backbone.Collection.extend({
  model: LeaderBoardUser,
  initialize: function() {
    this.api = 'users';
    this.days = 31;
  },
  url: function() { return '/api/leaderboard/'+ this.api +'/'+ this.days +'/?format=json'; },
});

/*
 * Views
 */
LeaderBoardUserView = Backbone.Marionette.View.extend({
  template: _.template('<p><% if(hover){ %><%= score %><% } else {%><%= name %><% } %></p>'),
  tagName: 'a',
  className: 'list-group-item',

  /* These events are only triggered when over
   * a span in the paragraph */
  events : {
    'mouseover': function() {
      var model = this.model;
      model.set('hover', true);
      setTimeout(function(){ model.set('hover', false); }, 5000); },
  },

  /* Setup event listeners for word spans */
  initialize : function(options) {
    this.listenTo(this.model, 'change:hover', this.render);
  },

  /* Triggers the proper class assignment
   * when the word <span> is redrawn */
  numberWithCommas: function(x) {
      return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",");
  },

  onRender : function() {
    if(this.model.get('user')) {
      this.$el.attr('href', '/u/'+this.model.get('name')+'/');
    } else {
      this.$el.attr('href', '/team/'+this.model.get('name')+'/');
    }

    if(this.model.get('hover')) {
      var $el = this.$el.find('p');
      $el.text(this.numberWithCommas($el.text()));
    }
  },
});


LeaderBoardUserListView = Backbone.Marionette.CollectionView.extend({
  childView  : LeaderBoardUserView,
  childViewEventPrefix: 'leaderboard:user',
});


LeaderBoardView = Backbone.Marionette.View.extend({
  template: '#dashboard-leaderboard-template',
  className: 'col-xs-2',

  regions: {
    'list': 'ol.list-group'
  },

  initialize: function() {
    var settings = new LeaderBoardSettings();
    this.collection = new LeaderBoardUserList();
    // var self = this;
    // this.$('h2').on('click', function() { self.model.toggle_api(); });
    // this.$('h4').on('click', function() { self.model.toggle_time_range(); });
    // this.listenTo(this.model, 'change:time_range', function() {
    //   var time_obj = this.model.time_key[this.model.get('time_range')];
    //   this.model.set('text', time_obj.text);
    //   this.model.set('days', time_obj.days);
    //   this.collection.days = time_obj.days;
    //   this.collection.fetch();
    // });
    // this.listenTo(this.model, 'change:api', function() {
    //   this.collection.api = this.model.get('api');
    //   this.collection.fetch();
    //   this.render();
    // });
    // this.listenTo(this.model, 'change:text', function() {
    //   this.render();
    // });
    this.collection.fetch();
  },

  onRender: function() {
    this.showChildView('list', new LeaderBoardUserListView({'collection': this.collection}) );
  }

});
