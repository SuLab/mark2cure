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
  template: '#relation-list-item-template',
  tagName: 'a',
  className: 'relation-item-container',

  onRender : function() {
    this.$el.attr('href', '/task/relation/'+ this.model.get('id') +'/')

    if(this.model.get('user').completed) {
      this.$el.addClass('disabled').attr('disabled', true).attr('href', '#');
    }
  },

});

DocumentCompositeView = Backbone.Marionette.CompositeView.extend({
  template: _.template('<div id="relation-task-list"></div>'),
  childView  : DocumentView,
  childViewContainer: "#relation-task-list",

  initialize : function(options) {
    this.collection.fetch();
  },

  onRender : function() {
  },

});

DocumentRelationBoard = new Backbone.Marionette.Application();
