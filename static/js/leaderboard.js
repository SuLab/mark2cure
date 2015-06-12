/*
 *  Models & Collections
 */
User = Backbone.RelationalModel.extend({
  defaults: {hover: false}
});

UserList = Backbone.Collection.extend({
  /* Common utils to perform on an array of Word
   * models for house keeping and search */
  model: User,
  url: '/api/leaderboard/users/?format=json',
});


/*
 * Views
 */
UserView = Backbone.Marionette.ItemView.extend({
  template: _.template('<p><a href="/u/<%- user.username %>" class="text-muted"><% if(hover){ %><%= rating_score %><% } else {%><%= user.username %><% } %></a></p>'),
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

UserCollectionView = Backbone.Marionette.CollectionView.extend({
  childView  : UserView,
  tagName   : 'ol',
  className : 'list-unstyled list-group',
  events : {
  },

  onRender : function() {
  },

});

LeaderBoard = new Backbone.Marionette.Application();
