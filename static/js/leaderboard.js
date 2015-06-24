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

User = Backbone.RelationalModel.extend({
  defaults: {'user': false, 'hover': false}
});

UserList = Backbone.Collection.extend({
  model: User,
  initialize: function() {
    this.api = 'users';
    this.days = 31;
  },
  url: function() { return '/api/leaderboard/'+ this.api +'/'+ this.days +'/?format=json'; },
});

/*
 * Views
 */
UserView = Backbone.Marionette.ItemView.extend({
  template: _.template('<p><a <% if(user){ %>href="/u/<%- name %>/"<% } else { %>href="/team/<%- name %>/"<% } %> class="text-muted"><% if(hover){ %><%= score %><% } else {%><%= name %><% } %></a></p>'),
  tagName: 'li',
  className: 'list-group-item',

  /* These events are only triggered when over
   * a span in the paragraph */
  events : {
    'mouseover': function() {
      var model = this.model;
      model.set('hover', true);
      setTimeout(function(){ model.set('hover', false); }, 1000); },
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
    if(this.model.get('hover')) {
      var $el = this.$el.find('a');
      $el.text(this.numberWithCommas($el.text()));
    }
  },

});

UserCompositeView = Backbone.Marionette.CompositeView.extend({
  template: _.template('<h2 class="text-center">Top <%- api %> <strong class="font-purple">></strong></h2><h4 class="text-center"><%- text %> <strong class="font-purple">></strong></h4><ol class="list-unstyled list-group"></ol>'),
  childView  : UserView,
  childViewContainer: "ol",

  initialize : function(options) {
    this.listenTo(this.model, 'change:time_range', function() {
      var time_obj = this.model.time_key[this.model.get('time_range')];
      this.model.set('text', time_obj.text);
      this.model.set('days', time_obj.days);
      this.collection.days = time_obj.days;
      this.collection.fetch();
    });

    this.listenTo(this.model, 'change:api', function() {
      this.collection.api = this.model.get('api');
      this.collection.fetch();
      this.render();
    });

    this.listenTo(this.model, 'change:text', function() {
      this.render();
    });

    this.collection.fetch();
  },

  onRender : function() {
    var self = this;
    this.$('h2').on('click', function() { self.model.toggle_api(); });
    this.$('h4').on('click', function() { self.model.toggle_time_range(); });
  },

});

LeaderBoard = new Backbone.Marionette.Application();
