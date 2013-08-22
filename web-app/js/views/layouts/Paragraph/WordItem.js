define(['marionette', 'templates', 'vent',
        //-- ETC
        'moment', 'd3'],
    function (Marionette, templates, vent) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.paragraph.worditem,
    tagName : "span",

    className : function() {
      if( this.model.get('parentDocument').get('complete') && 
          _.contains(this.options.ann_range, this.model.get('start')) ) { return 'selectedCompare' }
    },

    ui : {
      word : 'span.word'
    },

    events : {
      'click'     : 'clickWord',
      'mouseover' : 'hoveredWord',
      'mousemove' : 'hoveredWord',
      'mouseout'  : 'leaveWord',
    },

    initialize : function(options) {
      this.bindTo(this.model, 'change', this.render);
      // console.log('WORD INIT :: ', this);
    },

    onRender : function() {
      var self = this;
      this.renderingClassSetting('selected');
      this.renderingClassSetting('latest');

      if( this.model.get('parentDocument').get('complete') ) {
        var popularity = this.model.get('parentDocument').get('popularity')[ this.model.get('position') ];
        this.$el.css({'backgroundColor': this.options.color_scale(popularity)});
      }

    },

    //
    //-- Event actions
    //
    clickWord : function(evt) {
      evt.preventDefault();
      var self = this,
          model = self.model,
          collection = self.model.collection;
      if(model.get('parentDocument').get('complete')) { return false; }

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

    hoveredWord : function(evt) {
      evt.preventDefault();
      if(this.model.get('parentDocument').get('complete')) { return false; }
      this.ui.word.addClass('hover');
      this.updatePosition();
    },

    leaveWord : function(evt) {
      evt.preventDefault();
      this.ui.word.removeClass('hover');
    },

    //-- Utilities for view
    renderingClassSetting : function(attrCheck) {
      if( this.model.get(attrCheck) ) {
        this.ui.word.addClass(attrCheck);
      } else {
        this.ui.word.removeClass(attrCheck);
      }
    },

    updatePosition : function() {
      var position = this.ui.word.position();
      this.model.set('pos_x',   position.left );
      this.model.set('pos_y',   position.top );
      this.model.set('width',   this.ui.word.width() );
      this.model.set('height',  this.ui.word.height() );
    },

  });
});
