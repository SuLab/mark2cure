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
      this.listenTo(this.model, 'change:selected', this.render);
      this.listenTo(this.model, 'change:neighbor', this.render);

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

        self.createAnns(start_i, stop_i)
      } else {
        if( prexisting ) {
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
        console.log(ann.get('text'), " :: ", ann.get('length'), ann.get('start'));
      });
    },

    createAnns : function(start, stop) {
      var self = this,
          doc = this.model.get('parentDocument'),
          annotations = doc.get('annotations'),
          text = doc.get('text').substring(start, stop),
          //-- Get the "pure" text
          text = text.replace(/^[^a-z\d]*|[^a-z\d]*$/gi, '');
      var type = User.get('sel_mode');

      _.each(this.getIndicesOf(text, doc.get('text'), false), function(v) {
        annotations.create({
          kind      : 0,
          type      : type,
          text      : text,
          length    : text.length,
          start     : v,
          stop      : v+text.length
        });
      })

    },

    //-- Utilities for view
    getIndicesOf : function(searchStr, str, caseSensitive) {
      var startIndex = 0, searchStrLen = searchStr.length;
      var index, indices = [];
      if (!caseSensitive) {
          str = str.toLowerCase();
          searchStr = searchStr.toLowerCase();
      }
      while ((index = str.indexOf(searchStr, startIndex)) > -1) {
          indices.push(index);
          startIndex = index + searchStrLen;
      }
      return indices;
    },

    selectWordsOfAnnotations : function() {
      var ann_range =  this.model.get('parentDocument').get('annotations').map(function(m) {
        return {
          'start' : m.get('start'),
          'text' : m.get('text')
        }
      });

      this.model.collection.clear('selected');
      this.model.collection.each(function(word) {
        //- If the word is part of an annotation
        var found = _.filter(ann_range, function(ann) {
          return  word.get('text').length <= ann.text.length &&
                  word.get('start') >= ann.start &&
                  word.get('start') < (ann.start+ann.text.length);
        });

        if( found.length ) {
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
