define(['marionette', 'templates', 'vent',
        //-- ETC
        'moment'],
    function (Marionette, templates, vent) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.entity_tag.worditem,
    tagName : "span",

    events : {
      'click'     : 'clickWord',
    },

    initialize : function(options) {
      this.bindTo(this.model, 'change:selected', this.render);
    },

    onRender : function() {
      this.renderingClassSetting('selected');
    },

    //
    //-- Event actions
    //
    clickWord : function(evt) {
      evt.preventDefault();
      var self = this,
          model = self.model,
          collection = self.model.collection;

      if( evt.shiftKey ) {
        //-- Shift key was held, select all between the this and the "latest"
        var current_position = model.get('position'),
            latest_position = collection.findWhere({latest: true}).get('position'),
            selection = [current_position, latest_position],
            included_ids = _.range( _.min(selection), _.max(selection)+1 );
        _.each(included_ids, function(index) { collection.at(index).set('selected', true); });
      } else {
        model.set('selected', !this.model.get('selected'));
      }

      collection.each(function(word) { word.set('latest', false); });
      model.set('latest', true);
    },

    //-- Utilities for view
    renderingClassSetting : function(attrCheck) {
      if( this.model.get(attrCheck) ) {
        this.$el.addClass(attrCheck);
      } else {
        this.$el.removeClass(attrCheck);
      }
    }

  });
});
