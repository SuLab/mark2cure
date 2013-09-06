define(['marionette', 'templates', 'vent',
        //-- Models
        'models/User',
        //-- ETC
        'moment'],
    function (Marionette, templates, vent,
              User) {
  'use strict';

  return Marionette.ItemView.extend({
    template : templates.main.game.entity_tag.worditem,
    tagName : "span",

    events : {
      'mousedown'   : 'clickOrInitDrag',
      'mouseup'     : 'releaseDrag',
      'mouseover'   : 'hover',
      'mouseout'    : 'leaveHover'
    },

    initialize : function(options) {
      this.bindTo(this.model, 'change:selected', this.render);
      this.bindTo(this.model, 'change:neighbor', this.render);

      //-- (TODO) Why is this method so slow?
      // var doc = this.model.get('parentDocument');
      // this.listenTo(doc.get('annotations'), "add", this.selectWordsOfAnnotations, this);
      // this.listenTo(doc.get('annotations'), "remove", this.selectWordsOfAnnotations, this);
    },

    onRender : function() {
      this.renderingClassSetting('selected');
      this.renderingClassSetting('neighbor');
    },

    //
    //-- Event actions
    //
    hover : function(evt) {
      if(evt.which) {
        var last_model = this.model.collection.findWhere({latest: true}),
            sel = [last_model.get('start'), this.model.get('start')],
            range = [_.min(sel), _.max(sel)],
            highlight_list = this.model.collection.selectBetweenRange(range[0], range[1]+1);
        this.selectWordsOfAnnotations();

        _.each(highlight_list, function(word) { word.set('selected', true); });
      } else {
        if(this.model.get('selected')) {
          var anns = this.model.get('parentDocument').get('annotations'),
              ann_list = anns.findContaining(this.model.get('start'));
          if(ann_list.length > 0) {
            var types = _.uniq( _.map(ann_list, function(ann){ return ann.get('type'); }) );
            this.$el.tooltip('destroy');
            this.$el.tooltip({
              trigger: 'manual',
              title: types.join(', ')
            });
            this.$el.tooltip('show');
          }
        }
      }
    },

    leaveHover : function(evt) {
      this.$el.tooltip('destroy');
    },

    clickOrInitDrag : function() {
      this.model.collection.clear('latest');
      this.model.set('latest', true);
      this.model.set('selected', true);
    },

    releaseDrag : function(evt) {
      var self = this,
          last_model = this.model.collection.findWhere({latest: true}),
          doc = this.model.get('parentDocument'),
          dragged = last_model != this.model,
          annotations = doc.get('annotations'),
          ann_range =  annotations.getRange(),
          prexisting = _.contains(ann_range, this.model.get('start'));

      var type = User.get('sel_mode');

      if(dragged) {
        var sel = [last_model.get('start'), this.model.get('stop')],
            range = [_.min(sel), _.max(sel)],
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
          type      : type,
          text      : doc.get('text').substring(start_i, stop_i),
          length    : stop_i - start_i,
          start     : start_i,
          stop      : stop_i
        });

      } else {
        if( prexisting ) {
          //-- If the single annotation or range started on a prexisting annotation
          _.each(annotations.findContaining( this.model.get('start') ), function(ann) { ann.destroy(); })
        } else {
          //-- If the single annotation or range started on a prexisting annotation
          annotations.create({
            kind      : 0,
            type      : type,
            text      : self.model.get('text'),
            length    : self.model.get('length'),
            start     : self.model.get('start'),
            stop      : self.model.get('stop')
          });
        }
      }

      this.selectWordsOfAnnotations();
      this.selectNeighborsOfAnnotations();

      // console.log('/ / / / / / / / / / / /   ::   ', range);
      // _.each(annotations.models, function(ann) {
      //   console.log(ann.get('text'), " :: ", ann.get('length'))
      // });
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

    selectNeighborsOfAnnotations : function() {
      this.model.collection.clear('neighbor');
      var anns = this.model.get('parentDocument').get('annotations');
      this.model.collection.each(function(word, word_idx) {
        //-- Is the person to your right selected?
        var left_neighbor = word.collection.at(   word_idx - 1),
            right_neighbor = word.collection.at(  word_idx + 1);

        if(right_neighbor && !right_neighbor.get('selected') && word.get('selected')) {
          word.set('neighbor', true);
        }

        if(anns.exactMatch(word).length > 0 && left_neighbor) {
          if(left_neighbor.get('selected')) {
            left_neighbor.set('neighbor', true);
          }
          word.set('neighbor', true);
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
