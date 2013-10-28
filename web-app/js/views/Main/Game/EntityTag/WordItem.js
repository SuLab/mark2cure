define(['marionette', 'templates', 'vent',
        'models/User', 'moment'],
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
      'mouseout'    : function() { this.$el.tooltip('destroy'); }
    },

    initialize : function(options) {
      this.listenTo(this.model, 'change:selected', this.render);
      this.listenTo(this.model, 'change:neighbor', this.render);
      options['auto_select_all'] = true;
    },

    onRender : function() {
      this.renderingClassSetting('selected');
      this.renderingClassSetting('neighbor');
    },

    //
    //-- Event actions
    //
    hover : function(evt) {
      //-- If you're dragging with the m.25use down to make a large selection
      console.log(evt);
      if(evt.which) {

        var last_model = this.model.collection.findWhere({latest: true}),
            sel = [last_model.get('start'), this.model.get('start')],
            range = [_.min(sel), _.max(sel)],
            highlight_list = this.model.collection.selectBetweenRange(range[0], range[1]+1);
        this.selectWordsOfAnnotations();
        _.each(highlight_list, function(word) { word.set('selected', true); });

      } else {

        //-- If you hover over a word that is highlighted
        if(this.model.get('selected')) {
          var anns = this.model.get('parentDocument').get('annotations'),
              ann_list = anns.findContaining(this.model.get('start'));
          if(ann_list.length > 0) {
            //-- What annotation types does that particular word have? Show them w/ Bootstrap tooltip
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

    clickOrInitDrag : function() {
      //-- onmousedown we just set the word to be the latest so that we can refernce it later
      //-- whent he user releases after staying put or moving around
      this.model.collection.clear('latest');
      this.model.set({'latest': true, 'selected': true});
    },

    releaseDrag : function(evt) {
      //-- onmouseup from the user
      console.log('releaseDrag');
      var self = this,
          last_model = this.model.collection.findWhere({latest: true}),
          annotations = this.model.get('parentDocument').get('annotations'),
          ann_range =  annotations.getRange();

      var type = User.get('sel_mode');

      if( last_model != this.model ) {
        //
        //-- If the user just finished making a drag selection
        //
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

        self.createAnns(start_i, stop_i)
      } else {
        //
        //-- If it was a single click
        //
        if( _.contains(ann_range, this.model.get('start')) ) {
          //-- If the single annotation or range started on a prexisting annotation
          _.each(annotations.findContaining( this.model.get('start') ), function(ann) { ann.destroy(); })
        } else {
          //-- If the single annotation or range started on a prexisting annotation
          self.createAnns(self.model.get('start'), self.model.get('stop')+1);
        }
      }

      this.selectWordsOfAnnotations();
      this.selectNeighborsOfAnnotations();

      console.log('/ / / / / / / / / / / /');
      _.each(annotations.models, function(ann) {
        console.log(ann.get('text'), " || ", ann.get('start'), ann.get('length'), ann.get('stop'));
      });
    },

    createAnns : function(start, stop) {
      var self = this,
          type = User.get('sel_mode'),
          doc = this.model.get('parentDocument'),
          annotations = doc.get('annotations'),
          text = doc.get('text').substring(start, stop);

      if(this.options.auto_select_all) {
        //-- Get the "pure" text
        text = this.clean(text);
        _.each(this.getIndicesOf(text, doc.get('text'), false), function(v) {
          annotations.create({
            kind      : 0,
            type      : type,
            text      : text,
            length    : text.length,
            start     : v,
            stop      : v+text.length
          });
        });
      } else {
        annotations.create({
          kind      : 0,
          type      : type,
          text      : text,
          length    : text.length,
          start     : start,
          stop      : stop
        });
      }
    },

    //-- Utilities for view
    getIndicesOf : function(needle, haystack, caseSensitive) {
      var startIndex = 0,
          needleLen = needle.length,
          index,
          indices = [];

      if (!caseSensitive) {
          haystack = haystack.toLowerCase();
          needle = needle.toLowerCase();
      }

      while ((index = haystack.indexOf(needle, startIndex)) > -1) {
          var sliced = this.clean( haystack.substring(index - 1, index + needleLen + 1) );

          console.log(sliced, needle);
          if(sliced === needle) { indices.push(index); }
          startIndex = index + needleLen;
      }

      // console.log(indices);
      return indices;
    },

    clean : function(text) {
      return _.str.clean(text).replace(/^[^a-z\d]*|[^a-z\d]*$/gi, '');
    },

    selectWordsOfAnnotations : function() {
      var self = this,
          offset = 0,
          clean_word;
      var ann_range =  this.model.get('parentDocument').get('annotations').map(function(m) {
          return { 'start' : m.get('start'), 'stop' : m.get('stop') }
        });

      //-- Iterate over the words and see if they are contained within any of the documents annotations
      this.model.collection.clear('selected');
      this.model.collection.each(function(word) {
        // (TODO) Cache this!
        clean_word = self.clean( word.get('text') );
        //-- Get the offset so we know that finding the annotation in the word will work!
        if(word.get('text') !== clean_word) { offset = word.get('text').indexOf(clean_word); }

        //- If the word is within annotation
        var found = _.filter(ann_range, function(ann) { return  word.get('start')+offset >= ann.start && word.get('stop') <= ann.stop; });
        // console.log('Found :: ', found, ' :: ', word.attributes);
        if( found.length ) { word.set('selected', true); }
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
