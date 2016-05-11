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
  template: '#template',
  tagName: 'a',

  /*template: _.template('<p><%= title %></p>'),
  tagName: 'a',
  className: 'list-group-item',
  */

  onRender : function() {
    this.$el.attr('href', '/task/relation/'+ this.model.get('id') +'/')

    if(this.model.get('user').completed) {
      this.$el.addClass('disabled').attr('disabled', true).attr('href', '#');
    }
  },

});

DocumentCompositeView = Backbone.Marionette.CompositeView.extend({
  template: _.template('<div class="relations-container-inside"><div class="relations-title">Document Relation Tasks</div><div id="relation-list"></div></div>'),
  childView  : DocumentView,
  childViewContainer: "#relation-list",


  initialize : function(options) {
    this.collection.fetch();
  },

  onRender : function() {
  },

});

DocumentRelationBoard = new Backbone.Marionette.Application();
