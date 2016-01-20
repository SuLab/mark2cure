Document = Backbone.RelationalModel.extend({
  defaults: {}
});

DocumentList = Backbone.Collection.extend({
  model: Document,
  url: function() { return '/api/relationships/?format=json'; },
});

/*
 * Views
 */
DocumentView = Backbone.Marionette.ItemView.extend({
  template: _.template('<p><%= title %></p>'),
  tagName: 'a',
  className: 'list-group-item',

  onRender : function() {
    this.$el.attr('href', '/task/relation/'+ this.model.get('id') +'/')

    if(this.model.get('user').completed) {
      this.$el.addClass('disabled').attr('disabled', true).attr('href', '#');
    }
  },

});

DocumentCompositeView = Backbone.Marionette.CompositeView.extend({
  template: _.template('<ol class="list-unstyled list-group"></ol>'),
  childView  : DocumentView,
  childViewContainer: "ol",

  initialize : function(options) {
    this.collection.fetch();
  },

  onRender : function() {
  },

});

DocumentRelationBoard = new Backbone.Marionette.Application();
