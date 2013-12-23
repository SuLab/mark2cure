YPet = new Backbone.Marionette.Application();

//
//-- Models & Collections
//

Word = Backbone.RelationalModel.extend({
  defaults: {
    text      : '',
    length    : null,
    start     : null,
    stop      : null,
    latest    : false,
    selected  : false,
    neighbor  : false,
  }
});

Annotation = Backbone.RelationalModel.extend({
  defaults: {
    kind    : "e",
    text    : '',
    length  : null,
    start   : null,
    stop    : null,
  },
  sync : function() { return false; }
});

WordList = Backbone.Collection.extend({
  model   : Word,
  url     : '/api/v1/words',

  clear : function(attr) {
    return this.each(function(word) { word.set(attr, false); });
  },

  selectBetweenRange : function(start, stop) {
    return this.filter(function(word){
      return word.get('start') >= start && word.get('start') <= stop;
    });
  }
});

AnnotationList = Backbone.Collection.extend({
  model   : Annotation,
  url     : '/api/v1/annotations',

  getRange : function() {
    //-- Returns back an array of indexes that are occupied by annotations within the text
    // example, 6 letter annotation: [14, 15, 16, 17, 18, 19]
    var range = []
    this.each(function(annotation) {
      range.push( _.range(annotation.get('start'), annotation.get('stop')+1)  );
    })
    return _.uniq( _.flatten(range) );
  },

  findContaining : function(index) {
    return this.filter(function(annotation) {
      return index >= annotation.get('start') && index <= annotation.get('stop');
    });
  },

  exactMatch : function(word) {
    return this.filter(function(annotation) {
      return  word.get('start') == annotation.get('start') &&
              word.get('stop') == annotation.get('stop');
    });
  },

  add : function(ann) {
    //-- Prevent duplicate annotations from being submitted
    var isDupe = this.any(function(_ann) {
        return  _ann.get('text') === ann.get('text') &&
                _ann.get('start') === ann.get('start');
    });
    if (isDupe) { return false; }
    Backbone.Collection.prototype.add.call(this, ann);
  }
});

Paragraph = Backbone.RelationalModel.extend({
  defaults: {
    text : '',
  },

  relations: [{
    type: 'HasMany',
    key: 'annotations',

    relatedModel: Annotation,
    collectionType: AnnotationList,

    reverseRelation : {
      includeInJSON: true,
    }
  }, {
    type: 'HasMany',
    key: 'words',

    relatedModel: Word,
    collectionType: WordList,

    reverseRelation : {
      key : 'parentDocument',
      includeInJSON: false,
    }
  }],

  parseText : function() {
    var step = 0,
        word_obj,
        words = _.map(_.string.words( this.get('text') ), function(word) {
          word_obj = {
            'text'      : word,
            'length'    : word.length,
            'start'     : step,
            'stop'      : step + word.length
          }
          step = step + word.length + 1;
          return word_obj;
        });

      //-- Remove any words if they previously existed and add the new ones
      _.each(this.get('words'), function(word) { word.destroy(); });
      this.get('words').add(words);
  }

});

//
//-- Views
//

WordView = Backbone.Marionette.ItemView.extend({
  template : '#word-template',
  tagName : 'span',

  className : function() {
    var classArr = [],
        self = this;
    _.each(['selected', 'neighbor', 'latest'], function(v) {
      if(self.model.get(v)) { classArr.push(v) }
    });
    return classArr.join(' ');
  },

  events : {
    'mousedown'   : 'clickOrInitDrag',
    'mouseup'     : 'releaseDrag',
    'mouseover'   : 'hover',
  },

  initialize : function(options) {
    this.listenTo(this.model, 'change:selected', this.render);
    this.listenTo(this.model, 'change:neighbor', this.render);
    options['auto_select_all'] = false;
    options['firefox'] = navigator.userAgent.toLowerCase().indexOf("firefox") > -1;
  },

  onRender : function() {
    this.$el.attr('class', _.result(this, 'className'));
  },

  //
  //-- Event actions
  //
  clickOrInitDrag : function() {
    //-- onmousedown we just set the word to be the latest so that we can refernce it later
    //-- whent he user releases after staying put or moving around
    this.model.collection.clear('latest');
    this.model.set({'latest': true, 'selected': true});
  },

  hover : function(evt) {
    var dragging = this.options.firefox ? 0 : evt.which;
    //-- If you're dragging with the mouse down to make a large selection
    if( dragging ) {
      var last_model = this.model.collection.findWhere({latest: true}) || this.model.collection.at(0),
          sel = [last_model.get('start'), this.model.get('start')],
          range = [_.min(sel), _.max(sel)],
          highlight_list = this.model.collection.selectBetweenRange(range[0], range[1]+1);
      this.selectWordsOfAnnotations();
      _.each(highlight_list, function(word) { word.set('selected', true); });
    }
  },

  releaseDrag : function(evt) {
    //-- onmouseup from the user
    var last_model = this.model.collection.findWhere({latest: true})  || this.model.collection.at(0),
        annotations = this.model.get('parentDocument').get('annotations');

    //-- If the user just finished making a drag selection
    if( last_model != this.model ) {
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

      this.createAnns(start_i, stop_i)

    //-- If it was a single click
    } else {
      //-- If click was on preexisting annotation to delete it
      if( _.contains(annotations.getRange(), this.model.get('start')) ) {
        _.each(annotations.findContaining( this.model.get('start') ), function(ann) { ann.destroy(); })

      //-- If a new single annotation or range
      } else {
        this.createAnns(this.model.get('start'), this.model.get('stop')+1);
      }
    }

    this.selectWordsOfAnnotations();

    // console.log('/ / / / / / / / / / / /');
    // _.each(annotations.models, function(ann) {
      // console.log(ann.get('text'), " || ", ann.get('start'), ann.get('stop'));
    // });
  },

  createAnns : function(start, stop) {
    var doc = this.model.get('parentDocument'),
        annotations = doc.get('annotations'),
        text = doc.get('text').substring(start, stop),
        annotation = this.sanitizeAnnotation(text, start);

    // if(this.options.auto_select_all) {
      // _.each(this.getIndicesOf(text, doc.get('text'), false), function(v) {
      //   annotations.create({
      //     text      : text,
      //     length    : text.length,
      //     start     : v,
      //     stop      : v+text.length
      //   });
      // });
    // } else {
      annotations.create({
        text      : annotation.text,
        length    : annotation.text.length,
        start     : annotation.start,
        stop      : annotation.start + annotation.text.length
      });
    // }
  },

  //-- Utilities for view
  // getIndicesOf : function(needle, haystack, caseSensitive) {
  //   var startIndex = 0,
  //       needleLen = needle.length,
  //       index,
  //       indices = [];

  //   if (!caseSensitive) {
  //       haystack = haystack.toLowerCase();
  //       needle = needle.toLowerCase();
  //   }

  //   while ((index = haystack.indexOf(needle, startIndex)) > -1) {
  //       if(this.clean( haystack.substring(index - 1, index + needleLen + 1) ) === needle) { indices.push(index); }
  //       startIndex = index + needleLen;
  //   }
  //   return indices;
  // },

  sanitizeAnnotation : function(full_str, start) {
    //-- Return the cleaned string and the (potentially) new start position
    var str = _.str.clean(full_str).replace(/^[^a-z\d]*|[^a-z\d]*$/gi, '');
    // console.log(full_str, str);
    return {'text':str, 'start': start+full_str.indexOf(str)};
  },

  selectWordsOfAnnotations : function() {
    //-- Iterate over the words and see if they are contained within any of the documents annotations
    var self = this,
        annotations = this.model.get('parentDocument').get('annotations'),
        word_index_average = 0;

    //-- Before we select the words to highlight, remove all the preexisting ones
    this.model.collection.clear('selected');
    this.model.collection.clear('neighbor');

    _.each(annotations.models, function(ann) {
      var words = self.model.collection.filter(function(word) {
        word_index_average = (word.get('start') + word.get('stop')) / 2;
        // console.log(word.get('text'), " // ", word_index_average, word.get('start'), word.get('stop'),  " // ", ann.get('start'), ann.get('stop') );
        return ann.get('start') < word_index_average && word_index_average < ann.get('stop');
      });

      if(words.length) {
        _.each(words, function(word) { word.set('selected', true); })
        _.last(words).set('neighbor', true);
      }

    });
  }

});

WordCollectionView = Backbone.Marionette.CollectionView.extend({
  itemView: WordView,
  tagName : 'p',
  className : 'paragraph',
});
