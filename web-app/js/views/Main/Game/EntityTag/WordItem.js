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
        var self = this,
            last_model = this.model.collection.findWhere({latest: true}),
            doc = this.model.get('parentDocument'),
            dragged = last_model != this.model;

        if(dragged) {
          //-- (TODO) Account for reverse dragging
          var start_i = last_model.get('start'),
              stop_i = this.model.get('stop');
          doc.get('annotations').create({
            kind      : 0,
            type      : 'disease',
            text      : doc.get('text').substring(start_i, stop_i),
            length    : stop_i - start_i,
            start     : start_i,
            stop      : stop_i
          });

        } else {
          doc.get('annotations').create({
            kind      : 0,
            type      : 'disease',
            text      : self.model.get('text'),
            length    : self.model.get('length'),
            start     : self.model.get('start'),
            stop      : self.model.get('stop')
          });
        }

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
