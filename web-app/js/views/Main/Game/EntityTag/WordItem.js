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
      this.bindTo(this.model.get('parentDocument').get('annotations'), 'change', this.selectWordsOfAnnotations);
    },

    onRender : function() {
      this.renderingClassSetting('selected');
    },

    //
    //-- Event actions
    //
    clickOrInitDrag : function() {
      this.model.collection.clear('latest');
      this.model.set('latest', true);
    },

    releaseDrag : function(evt) {
      var self = this,
          last_model = this.model.collection.findWhere({latest: true}),
          doc = this.model.get('parentDocument'),
          dragged = last_model != this.model,
          annotations = doc.get('annotations'),
          ann_range =  annotations.getRange();

        if(dragged) {
          var sel = [last_model.get('start'), this.model.get('stop')],
              range = Array.prototype.slice.call(sel).sort(),
              start_i = range[0],
              stop_i = range[1]+1;

          //-- (TODO) http://stackoverflow.com/questions/7837456/comparing-two-arrays-in-javascript
          if( String(sel) !== String(range) ) {
            //-- They dragged in reverse
            start_i = this.model.get('start');
            stop_i = last_model.get('stop')+1;
          }

          annotations.create({
            kind      : 0,
            type      : 'disease',
            text      : doc.get('text').substring(start_i, stop_i),
            length    : stop_i - start_i,
            start     : start_i,
            stop      : stop_i
          });

        } else {

          if( _.contains(ann_range, this.model.get('start')) ) {
            //-- If the single annotation or range started on a prexisting annotation
            _.each(annotations.findContaining( this.model.get('start') ), function(ann) { ann.destroy(); })
          } else {
            //-- If the single annotation or range started on a prexisting annotation
            annotations.create({
              kind      : 0,
              type      : 'disease',
              text      : self.model.get('text'),
              length    : self.model.get('length'),
              start     : self.model.get('start'),
              stop      : self.model.get('stop')
            });
          }

        }
      this.selectWordsOfAnnotations();
    },

    //-- Utilities for view
    selectWordsOfAnnotations : function() {
      var ann_range =  this.model.get('parentDocument').get('annotations').getRange();
      this.model.collection.clear('selected');
      this.model.collection.each(function(word) {
        //- If the word is part of an annotation
        if( _.contains(ann_range, word.get('start')) ) {
          word.set('selected', true);
        }
      });
    },

    renderingClassSetting : function(attrCheck) {
      if( this.model.get(attrCheck) ) {
        this.$el.addClass(attrCheck);
      } else {
        this.$el.removeClass(attrCheck);
      }
    }

  });
});
