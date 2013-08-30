define(['marionette', 'templates', 'vent',
        //-- ETC
        'moment'],
    function (Marionette, templates, vent) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.entity_tag.worditem,
    tagName : "span",

    events : {
      'mousedown'   : 'clickOrInitDrag',
      'mouseup'     : 'releaseDrag'
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
    clickOrInitDrag : function() {
      if(false) {
        //-- if this is part of a bigger, prexisting annotation?
      } else {
        this.model.set('selected', !this.model.get('selected'));
      }

      this.model.collection.each(function(word) { word.set('latest', false); });
      this.model.set('latest', true);
    },

    releaseDrag : function(evt) {
        var last_model = this.model.collection.findWhere({latest: true});
        console.log(last_model.get('text'), this.model.get('text'));
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
