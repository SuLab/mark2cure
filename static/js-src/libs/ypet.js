Date.now = Date.now || function() { return +new Date; };
(function (d) {
    d.each(["backgroundColor", "borderBottomColor", "borderLeftColor", "borderRightColor", "borderTopColor", "color", "outlineColor"], function (f, e) {
            d.fx.step[e] = function (g) {
                        if (!g.colorInit) {
                                        g.start = c(g.elem, e);
                            g.end = b(g.end);
                            g.colorInit = true
                                        }
                        g.elem.style[e] = "rgb(" + [Math.max(Math.min(parseInt((g.pos * (g.end[0] - g.start[0])) + g.start[0]), 255), 0), Math.max(Math.min(parseInt((g.pos * (g.end[1] - g.start[1])) + g.start[1]), 255), 0), Math.max(Math.min(parseInt((g.pos * (g.end[2] - g.start[2])) + g.start[2]), 255), 0)].join(",") + ")"
                    }
        });

    function b(f) {
            var e;
            if (f && f.constructor == Array && f.length == 3) {
                        return f
                    }
            if (e = /rgb\(\s*([0-9]{1,3})\s*,\s*([0-9]{1,3})\s*,\s*([0-9]{1,3})\s*\)/.exec(f)) {
                        return [parseInt(e[1]), parseInt(e[2]), parseInt(e[3])]
                    }
            if (e = /rgb\(\s*([0-9]+(?:\.[0-9]+)?)\%\s*,\s*([0-9]+(?:\.[0-9]+)?)\%\s*,\s*([0-9]+(?:\.[0-9]+)?)\%\s*\)/.exec(f)) {
                        return [parseFloat(e[1]) * 2.55, parseFloat(e[2]) * 2.55, parseFloat(e[3]) * 2.55]
                    }
            if (e = /#([a-fA-F0-9]{2})([a-fA-F0-9]{2})([a-fA-F0-9]{2})/.exec(f)) {
                        return [parseInt(e[1], 16), parseInt(e[2], 16), parseInt(e[3], 16)]
                    }
            if (e = /#([a-fA-F0-9])([a-fA-F0-9])([a-fA-F0-9])/.exec(f)) {
                        return [parseInt(e[1] + e[1], 16), parseInt(e[2] + e[2], 16), parseInt(e[3] + e[3], 16)]
                    }
            if (e = /rgba\(0, 0, 0, 0\)/.exec(f)) {
                        return a.transparent
                    }
            return a[d.trim(f).toLowerCase()]
        }
    function c(g, e) {
            var f;
            do {
            f = d.css(g, e);
            if (f != "" && f != "transparent" || d.nodeName(g, "body")) {
                            break
                        }
                        e = "backgroundColor"
        } while (g = g.parentNode);
                return b(f)
                    }
    var a = {
            aqua: [0, 255, 255],
            azure: [240, 255, 255],
            beige: [245, 245, 220],
            black: [0, 0, 0],
            blue: [0, 0, 255],
            brown: [165, 42, 42],
            cyan: [0, 255, 255],
            darkblue: [0, 0, 139],
            darkcyan: [0, 139, 139],
            darkgrey: [169, 169, 169],
            darkgreen: [0, 100, 0],
            darkkhaki: [189, 183, 107],
            darkmagenta: [139, 0, 139],
            darkolivegreen: [85, 107, 47],
            darkorange: [255, 140, 0],
            darkorchid: [153, 50, 204],
            darkred: [139, 0, 0],
            darksalmon: [233, 150, 122],
            darkviolet: [148, 0, 211],
            fuchsia: [255, 0, 255],
            gold: [255, 215, 0],
            green: [0, 128, 0],
            indigo: [75, 0, 130],
            khaki: [240, 230, 140],
            lightblue: [173, 216, 230],
            lightcyan: [224, 255, 255],
            lightgreen: [144, 238, 144],
            lightgrey: [211, 211, 211],
            lightpink: [255, 182, 193],
            lightyellow: [255, 255, 224],
            lime: [0, 255, 0],
            magenta: [255, 0, 255],
            maroon: [128, 0, 0],
            navy: [0, 0, 128],
            olive: [128, 128, 0],
            orange: [255, 165, 0],
            pink: [255, 192, 203],
            purple: [128, 0, 128],
            violet: [128, 0, 128],
            red: [255, 0, 0],
            silver: [192, 192, 192],
            white: [255, 255, 255],
            yellow: [255, 255, 0],
            transparent: [255, 255, 255]
        }
})(jQuery);


/*
 *  Models & Collections
 */

NERAnnotationTypeList = Backbone.Collection.extend({
  /* Very simple collection to store the type of
   * Annotations that the application allows
   * for paragraphs */
  model: Backbone.Model.extend({}),

  url: function() { return false; },
  sync: function () { return false; },
});

NERAnnotationTypes = new NERAnnotationTypeList([
  {name: 'Disease', color: '#d1f3ff'},
  {name: 'Gene', color: '#B1FFA8'},
  {name: 'Drug', color: '#ffd1dc'}
]);

NERWord = Backbone.RelationalModel.extend({
  /* A Word model represents each tokenized word present
   * in the paragraph YPet is attached to. */
  defaults: {
    text: '',
    start: null,
    latest: null,
    neighbor: false,
  },

  url: function() { return false; },
  sync: function () { return false; }
});

NERWordList = Backbone.Collection.extend({
  /* Common utils to perform on an array of Word
   * models for house keeping and search */
  model: NERWord,

  url: function() { return false; },
  sync: function () { return false; }
});

NERAnnotation = Backbone.RelationalModel.extend({
  /* Each annotation in the paragraph. An Annotation
   * is composed of an array of Words in order to determine
   * the full text and position which are not
   * explicity set */
  defaults: {
    /* An annotation doesn't exist when removed so
     * we can start them all off at 0 and not need to
     * mix in a null type */
    type_id: 0,
    text: '',
    start: null,
    self: true, // To determine if they're the author's or an opponent's
  },

  url: function() { return false; },
  sync: function () { return false; },

  relations: [{
    type: 'HasMany',
    key: 'words',

    relatedModel: NERWord,
    collectionType: NERWordList
  }],

  validate: function(attrs) {
    if(attrs.text == null) {
      return 'We need text!'
    }
  },

  initialize: function() {
    var self = this;
    /* Santize BioC passage or other sources of
     * Annotation information into the Backbone model
     */
    this.set('id', +this.get('@id'));
    this.unset('@id');
    if(this.get('infon')) {
      _.each(this.get('infon'), function(i) {
        if(!isNaN(i['#text'])) {
          self.set(i['@key'], +i['#text'])
        } else {
          self.set(i['@key'], i['#text'])
        }
      });
      this.unset('infon');
    }
    if(this.get('location')) {
      this.set('start', +this.get('location')['@offset']);
      if(!this.get('text')){
        return false;
      }
      if(this.get('text').length != +this.get('location')['@length']) {
        Raven.captureMessage('Text length and reported @length failure', {extra: {'ann_db_pk': this.get('uid')}});
      }
      this.unset('location');
    }

    /* Find the word instances needed from the
     * instances from the paragraph */
    var selected_words = this.get('paragraph').get('words').filter(function(word) {
      return word.get('start') >= self.get('start') && word.get('start') < self.get('start')+self.get('text').length;
    });
    /*
      var ann_start = +annotation.location['@offset'] - passage_offset;
      var ann_length = +annotation.location['@length'];
      var ann_type_id = +_.find(annotation.infon, function(o){return o['@key']=='type_id';})['#text'];

      var start_match = false;
      var selected = words.filter(function(word) {
        // The Annotation found a word which matches start position exactly
        var starts = word.get('start') == ann_start;
        if (starts) { start_match = true; }
        return starts || ( word.get('start') > ann_start && word.get('start') < ann_start+ann_length );
      });

      try {
        ...
        var words_match = selected.length == _.str.words(annotation.text).length;
        if(words_match==false && start_match==false) {
          Raven.captureMessage('Imperfect Pubtator >> YPet Match', {extra: {
            'selected': selected,
            'annotation': annotation,
            'passage': passage
          }});
        }

      } catch(e) { Raven.captureException(e); }


    */

    this.get('words').each(function(word) { word.destroy(); });
    this.get('words').set(selected_words);

    var words_len = this.get('words').length;
    this.get('words').each(function(word, word_index) {
      if(word_index == words_len-1) { word.set('neighbor', true); }
      word.trigger('highlight');
    });
  },

  toggleType: function() {
    /* Removes (if only 1 Annotation type) or changes
     * the Annotation type when clicked after existing */
    if( this.get('type_id') == NERAnnotationTypes.length-1 || this.get('text') == "") {
      this.destroy();
    } else {
      this.set('type_id', this.get('type_id')+1 );
    }
  }
});


NERAnnotationList = Backbone.Collection.extend({
  /* Utils for the Paragraph Annotations lists
   * collectively */
  model: NERAnnotation,

  url: function() { return false; },
  sync: function () { return false; },

  sanitizeAnnotation : function(full_str, start) {
    /* Return the cleaned string and the (potentially) new start position */
    var str = _.str.clean(full_str).replace(/^[^a-z\d]*|[^a-z\d]*$/gi, '');
    return {'text':str, 'start': start+full_str.indexOf(str)};
  },

  initialize: function(options) {
    this.listenTo(this, 'add', function(annotation) {
      var ann = this.sanitizeAnnotation(annotation.get('words').pluck('text').join(' '), annotation.get('words').first().get('start'));
      annotation.set('text', ann.text);
      annotation.set('start', ann.start);
      this.drawAnnotation(annotation);
    });
    this.listenTo(this, 'change:type_id', function(annotation) {
      this.drawAnnotation(annotation);
    });
    this.listenTo(this, 'remove', function(annotation, collection) {
      /* Must iterate b/c annotation var "words" attribute is
       * empty at this point */
      collection.parentDocument.get('words').each(function(word) {
        word.trigger('highlight', {'color': '#fff'});
        word.set('neighbor', false);
      });

      collection.each(function(annotation) {
        collection.drawAnnotation(annotation);
      });
    });
  },

  drawAnnotation: function(annotation) {
    /* this = NERAnnotationList collection */
    var annotation_type = NERAnnotationTypes.at(annotation.get('type_id'));
    // var parent_document = this.parentDocument || this._parentDocument; // undefined

    /* Draw all the basic background or underlines */
    annotation.get('words').each(function(word, word_index) {
      // if(annotation.get('opponent')) {
        // word.trigger('underline', {'color': annotation_type.get('color')});
      // } else {
      word.trigger('highlight', {'color': annotation_type.get('color')});
      // }
    });

    if(annotation.get('opponent')) {
      var words = annotation.get('words')
      var author_annotations = parent_document.get('annotations');

      var anns = []
      author_annotations.each(function(main_ann) {
        if(main_ann.get('words').contains(words.first()) || main_ann.get('words').contains(words.last())) {
          anns.push(main_ann.cid);
        }
      });

      if(_.uniq(anns).length > 1) {
        /* 2 Different Parent Annotations */
        words.each(function(word) {
          if(word == words.last()) {
            word.trigger('underline-space', {'color': '#fff', 'last_word': true});
          } else {
            word.trigger('underline-space', {'color': annotation_type.get('color'), 'last_word': false});
          }
        });
      }
    }
  }

  /*
    var user_ids = this.model.get('annotations').pluck('user_id');
    //-- Why would these annotations contain multiple user_ids?
    if(user_ids.length >= 1) { console.log('throw error'); }
    var user_id = +user_ids[0];
  */

});


NERParagraph = Backbone.RelationalModel.extend({
  /* Foundational model for tracking everything going
   * on within a Paragraph like Words and Annotations */
  defaults: {
    text: '',
    offset: 0,
  },

  url: function() { return false; },
  sync: function () { return false; },

  relations: [{
    /* Many Words */
    type: 'HasMany',
    key: 'words',

    relatedModel: NERWord,
    collectionType: NERWordList,

    reverseRelation : {
      key : 'paragraph',
      includeInJSON: false,
    }
  }, {
    /* Many Annotations */
    type: 'HasMany',
    key: 'annotations',

    relatedModel: NERAnnotation,
    collectionType: NERAnnotationList,

    reverseRelation : {
      key : 'paragraph',
      includeInJSON: false,
    }
  }],

  initialize : function() {
    /* Extract (tokenize) the individual words */
    var self = this;
    var step = 0,
        space_padding,
        word_obj,
        text = this.get('text'),
        words = _.map(_.str.words( text ), function(word) {
          word_obj = {
            'text': word,
            'start': step
          }
          space_padding = (text.substring(step).match(/\s+/g) || [""])[0].length;
          step = step + word.length + space_padding;
          return word_obj;
        });
    this.set('words', new NERWordList(words));

    // Normalize BioC values
    this.set('offset', +this.get('offset'));
    if(this.get('infon')) {
      _.each(this.get('infon'), function(i) {
        if(!isNaN(i['#text'])) {
          self.set(i['@key'], +i['#text'])
        } else {
          self.set(i['@key'], i['#text'])
        }
      });
      this.unset('infon');
    }

    /* Set the Annotations correctly for this paragraph */
    _.each(this.get('annotation'), function(obj) { obj.paragraph = self; });
    _.each(this.get('annotations'), function(obj) { obj.paragraph = self; });
    if(this.get('annotation') && Array.isArray(this.get('annotation'))) {
      this.set('annotations', new NERAnnotationList(this.get('annotation')));
      this.unset('annotation');
    }
  },
});

NERParagraphList = Backbone.Collection.extend({
  model: NERParagraph,

  url: function() { return false; },
  sync: function () { return false; },
});

NERDocument = Backbone.RelationalModel.extend({
  /* A Document with all the children required
   * for NER tasks such as:
   * - Passage offset generation
   * - Passage word tokenization
   * - Passage user_annotations
   * - Passage opponent_annotations
   */
  defaults: {
    pk: 0,
    pmid: '',
    active: false, // If is the currently addressed Document in the Quest
    completed: false, // If has previously been completed by the user
  },

  relations: [{
    type: 'HasMany',
    key: 'passages',

    relatedModel: NERParagraph,
    collectionType: NERParagraphList
  }]
});

NERDocumentList = Backbone.Collection.extend({
  /* The Quest of Documents to Complete */
  model: NERDocument,

  url: function() { return false; },
  sync: function () { return false; },
})


/*
 * Views
 */

NERWordView = Backbone.Marionette.View.extend({
  /* View for all direct actions on a word
  * - Model = NERWord
  * - Collection = None
  */
  template: _.template('<% if(neighbor) { %><%= text %><% } else { %><%= text %> <% } %>'),
  tagName: 'span',

  modelEvents: {
    'change:neighbor': function() { this.render(); },
    'change:disabled': 'actOnChangeDisabled',
    'change:latest': 'actOnChangeLatest',
    'change:masked': 'actOnChangeMasked',

    'unclick': 'actOnUnClick',
    'highlight': 'actOnHighlight',
    'underline': 'actOnUnderline',
    'underline-space': 'actOnUnderlineSpace'
  },

  actOnChangeLatest: function(model, value, options) {
    if(this.model.get('latest')) {
      this.model.trigger('highlight', {'color': '#D1F3FF'});
    }
    if(options.force) {
      this.model.trigger('highlight', {'color': '#fff'});
    }
  },

  actOnChangeDisabled: function() {
    if(this.model.get('disabled')) {
      this.$el.css('cursor', 'not-allowed');
    }
  },

  actOnHighlight: function(options) {
    this.$el.css({'backgroundColor': options.color});
  },

  actOnUnClick: function() {
    var $el = this.$el;
    $el.animate({backgroundColor: '#fff'}, 750, function() {
      $el.trigger('mousedown').trigger('mouseup');
    });
  },

  actOnChangeMasked: function() {
    if(this.model.get('masked')) {
      this.$el.css({'color': '#000', 'cursor': 'default', 'opacity': '.5'});
    }
  },

  actOnUnderline: function() {
    var $container = this.$el.parent(),
        pos = this.$el.position(),
        split_end = this.$el.height() >= 30; /* (TODO) Compare to reference single height unit */

    var yaxis = pos.top + this.$el.height() + 2;
    var width = this.$el.width() + 1;

    if (split_end) {
      /* The first part of the word that wraps to the second line */
      var absolute_left = $container.find('span').first().position().left;
      var split_left = $prev.position().left + $prev.width();
      var $prev = this.$el.prev(),
          $next = this.$el.next();

      $container.append('<div class="underline" style=" \
        position: absolute; \
        height: 4px; \
        width: '+ (Math.abs( pos.left+width - split_left)) +'px; \
        top: '+ (pos.top+(this.$el.height()/2)-5)  +'px; \
        left: '+ split_left +'px; \
        background-color: '+ d3.rgb(options.color).darker(.5) +';"></div>');

      /* The reminder on the line below */
      /* (TODO) sometimes it'll split and there will be no next word */
      $container.append('<div class="underline" style=" \
        position: absolute; \
        height: 4px; \
        width: '+ ($next.position().left - absolute_left) +'px; \
        top: '+ yaxis +'px; \
        left: '+ absolute_left +'px; \
        background-color: '+ d3.rgb(options.color).darker(.5) +';"></div>');

    } else {
      $container.append('<div class="underline" style=" \
        position: absolute; \
        height: 4px; \
        width: '+ width +'px; \
        top: '+ yaxis +'px; \
        left: '+ pos.left +'px; \
        background-color: '+ d3.rgb(options.color).darker(.5) +';"></div>');
    }
  },

  actOnUnderlineSpace: function() {
    var $container = this.$el.parent(),
    pos = this.$el.position(),
    color = d3.rgb(options.color).darker(2);

    var yaxis = pos.top + this.$el.height() + 2;
    var width = this.$el.width();
    if(options.last_word) {
      width = width - 5;
      color = '#fff';
    }

    $container.append('<div class="underline-space" style=" \
      position: absolute; \
      height: 4px; \
      width: 5px; \
      top: '+ yaxis +'px; \
      left: '+ (pos.left+width) +'px; \
      background-color: '+ color +';"></div>');
  },

  /* Triggers the proper class assignment
   * when the word <span> is redrawn */
  onRender : function() {
    this.$el.css(this.model.get('neighbor') ?
      {'margin-right': '5px', 'padding-right': '0'} :
      {'margin-right': '0px', 'padding-right': '.25rem'});
  },

  /* These events are only triggered when over
   * a span in the paragraph */
  events : {
    'mousedown' : 'mousedown',
    'mouseover' : 'mouseover',
    'mouseup'   : 'mouseup',
  },

  mousedown : function(evt) {
    /* Upon first clicking down, flag that entry
     * word as the start of the annotation */
    evt.stopPropagation();
    if(this.model.get('disabled')) { return; };
    this.model.set({'latest': 1});
  },

  mouseover : function(evt) {
    /* When selecting the annotation is ongoing
     * This is likely during dragging over words
     * to make a span selection */
    evt.stopPropagation();
    if(this.model.get('disabled')) { return; };
    var word = this.model,
        words = word.collection;

    /* You're dragging if another word has a latest timestamp */
    if(_.compact(words.pluck('latest')).length) {
      if(_.isNull(word.get('latest'))) { word.set({'latest': Date.now()}); }

      /* If the hover doesn't proceed in ordered fashion
       * we need to "fill in the blanks" between the words */
      var current_word_idx = words.indexOf(word);
      var first_word_idx = words.indexOf( words.find(function(word) { return word.get('latest') == 1; }) );

      /* Select everything from the starting to the end without
       * updating the timestamp on the first_word */
      var starting_positions = first_word_idx <= current_word_idx ? [first_word_idx, current_word_idx+1] : [first_word_idx+1, current_word_idx];
      var selection_indexes = _.range(_.min(starting_positions), _.max(starting_positions));
      _.each(_.without(selection_indexes, first_word_idx), function(idx) { words.at(idx).set('latest', Date.now()); });

      /* If there are extra word selections up or downstream
       * from the current selection, remove those */
      var last_selection_indexes = _.map(words.reject(function(word) { return _.isNull(word.get('latest')); }), function(word) { return words.indexOf(word); });
      var remove_indexes = _.difference(last_selection_indexes, selection_indexes);

      var word, ann;
      _.each(remove_indexes, function(idx) {
        word = words.at(idx);
        word.set('latest', null, {'force': true});
        // (TODO) What is going on here?
        ann = word.get('parentAnnotation');
        if(ann) { ann.collection.drawAnnotation(ann); }
      });
    }
  },

  mouseup : function(evt) {
    /* When selecting the annotation has stopped
     * This could be after a click, drag, or any other
     * event suggesting the selection is over. */
    evt.stopPropagation();
    if(this.model.get('disabled')) { return; };
    var word = this.model,
        words = word.collection;

    var selected = words.filter(function(word) { return word.get('latest') });
    if(selected.length == 1 && word.get('parentAnnotation') ) {
      word.get('parentAnnotation').toggleType();
    } else {
      /* if selection includes an annotation, delete that one */
      _.each(selected, function(w) {
        if(w.get('parentAnnotation')) {
          w.get('parentAnnotation').destroy();
        }
      });
      /* Take the selected words as a whole and make a new Annotation with them */
      console.log(word);
      word.get('parentDocument').get('annotations').create({words: selected});
    };

    words.each(function(word) { word.set('latest', null); });
  }
});

NERWordsView = Backbone.Marionette.CollectionView.extend({
  /* List of individual words
   * this.model = None
   * this.collection = NERWordList
  */
  tagName: 'p',
  className: 'paragraph m-0',
  childView: NERWordView,
  childViewEventPrefix: 'word'
});

NERParagraphView = Backbone.Marionette.View.extend({
  /* Container item for the individual words
   * this.model = NERParagraph
   * this.collection = None
  */
  template: _.template('<div class="paragraph-content p-3"></div>'),
  className: function() {
    var enabled_status = this.getOption('enabled') ? 'ypet-enabled' : 'ypet-disabled';
    return ['paragraph-box', 'mb-2', enabled_status].join(' ');
  },

  regions: {
    'words': '.paragraph-content'
  },

  onRender : function() {
    if(!this.options.enabled) {
      /* If you're showing a partner's results, disallow highlighting */
      // this.$el.css({'color': '#000', 'cursor': 'default'});;
      this.model.get('words').each(function(w) {
        w.set('disabled', true);
      });
      this.$el.css('cursor', 'not-allowed');

      console.log('onRender of NER Para View', this.$el);
    }

    this.showChildView('words', new NERWordsView(
      {'collection': this.model.get('words'),
        'enabled': this.options.enabled
      }));

    var has_opponent =  _.contains(this.model.get('annotations').pluck('self'), false);
    if (has_opponent) {
      /* Warning: this allows paragraph specific disabiling */
      this.triggerMethod('view:only');
    }
  },

  onViewOnly: function() {
  },

  events: {
    // 'mousedown': 'startCapture',
    // 'mousemove': 'startHoverCapture',
    // 'mouseup': 'captureAnnotation',
    // 'mouseleave': 'captureAnnotation',
  },
});

NERParagraphsView = Backbone.Marionette.CollectionView.extend({
  /* Parent list for NERParagraphs
   * this.model = NERDocument
   * this.collection = NERParagraphList
   */
  childView: NERParagraphView,
  className: 'paragraphs',
  childViewEventPrefix: 'paragraph',

  childViewOptions: function() {
    var enabled = true;
    if(this.options.mode == 'er') { enabled = false; }
    return {
      'enabled': enabled
    }
  },
  initialize: function() {
    console.log('ner para init', this);
  }
});

NERLoadingView = Backbone.Marionette.View.extend({
  /* Initial HTML before a REExtractionList is available
  * - Model = None
  * - Collection = None
  */
  template: '<p>Text Loading...</p>'
});

NERNavigationView = Backbone.Marionette.View.extend({
  /* The Progress indicator view for all interactions on the Quest
   * this.model = None
   * this.collection = NERDocumentList (Quest Documents)
   */
  template: '#ypet-navigation-template',
  className: 'row',

  onRender: function() {
     // #quest-guide.row
     //  .col-xs-10.col-xs-offset-1.col-md-5.col-md-offset-1.col-lg-3.col-lg-offset-1
     //    ol.list-unstyled.list-inline
     //      - for doc in completed_doc_pks
     //        li.list-inline-item.completed &#8226;
     //
     //      - each val, index in uncompleted_docs
     //        - if index == 0
     //          li.list-inline-item.active &#8226;
     //        - else
     //          li.list-inline-item.uncompleted &#8226;
     //
     //  .col-xs-10.col-xs-offset-1.col-md-3.col-md-offset-2.col-lg-3.col-lg-offset-4
     //    p.text-xs-center Score: <span id='score'>#{user.userprofile.score}</span>
  }
});

NERFooterView = Backbone.Marionette.View.extend({
});

YPet = Backbone.Marionette.View.extend({
  /* The top level view for all interactions on text
   * this.model = active NERDocument being annotated
   * this.collection = NERDocumentList (Quest Documents)
   */
  template: '#ypet-template',
  className: 'row',

  regions: {
    'navigation': '#ypet-navigation',
    'text': '#ypet-text',
    'footer': '#ypet-footer'
  },

  initialize: function() {
    var self = this;

    if(!self.options.training && !self.collection) {
      $i.getJSON('/document/'+self.options.pmid+'/', function( data ) {
        /* The Annotation information has been returned from the server at this point
           it is now safe to start YPET */
        self.model = new NERDocument(data);
        self.render();
      });
    }
  },

  onRender: function() {
    console.log('YPet onRender', this);

    if(this.collection && this.model) {
      this.options['collection'] = this.collection;
      this.showChildView('text', new NERParagraphsView(this.options));

      if(this.options.mode == 'ner') {
        this.showChildView('navigation', new NERNavigationView(this.options));
        this.showChildView('footer', new RelationTextView(this.options));

      } else if (this.options.mode == 're') {
        // var concept_uids = [this.options.concepts['c1'].id, this.options.concepts['c2'].id];
      }
    } else {
      this.showChildView('text', new NERLoadingView());
    };


      // var passages;
      //
      // YPet.addInitializer(function(options) {
      //
      //   #<{(| Kick off the application with the document data
      //   * - Even if pubtator is down / no annotations available
      //   * - the passage text is still available to YPet
      //   |)}>#
      //   $.getJSON('/document/pubtator/{{document.document_id}}.json', function( data ) {
      //     passages = data.collection.document.passage;
      //     var regions = {};
      //
      //     _.each(passages, function(passage, passage_idx) {
      //       var passage_id = _.find(passage.infon, function(o){return o['@key']=='id';})['#text'];
      //       var p_body = '<div id="'+ passage_id +'" class="paragraph-box m-t-1"><p class="paragraph"></p></div></div>';
      //       $('.paragraphs').append(p_body);
      //       regions[''+passage_idx] = '#'+passage_id;
      //     });
      //     YPet.addRegions(regions);
      //
      //     #<{(| Add paragraph words and pubtator annotations if available |)}>#
      //     _.each(passages, function(passage, passage_idx) {
      //       var p = new Paragraph({'text': passage.text});
      //       YPet[''+passage_idx].show( new WordCollectionView({
      //         collection: p.get('words'),
      //         passage_json: passage,
      //         bioc_json: data
      //       }) );
      //     });
      //
      //     #<{(| Add event listener for recent annotation search helper |)}>#
      //     _.each(passages, function(passage, passage_idx) {
      //       YPet[''+passage_idx].currentView.collection.parentDocument.get('annotations').on('change', function(model, collection) {
      //         var model_json = model.toJSON();
      //         set_google_annotation(model_json.text);
      //       });
      //     });
      //
      //   });
      // });
      // YPet.start();
      //
      // function set_google_annotation(ann_text) {
      //   var url = 'https://www.google.com/search?q='+ann_text;
      //   $('#google_annotation a').attr('href', url);
      //   $('#google_annotation a small').text(_.str.truncate(ann_text, 36));
      // };
      //
      // function show_results() {
      //   #<{(| Show partner comparison |)}>#
      //   $.ajax({
      //     url:'/task/entity-recognition/{{task.pk}}/{{document.pk}}/results.json',
      //     dataType: 'json',
      //     cache: false,
      //     async: false,
      //     success: function(data) {
      //
      //       window.scrollTo(0,0);
      //       var points = _.find(data.collection.infon, function(o){return o['@key']=='points';})['#text'];
      //       var partner = _.find(data.collection.infon, function(o){return o['@key']=='partner';})['#text'];
      //       var partner_level = _.find(data.collection.infon, function(o){return o['@key']=='partner_level';})['#text'];
      //       var flatter = _.find(data.collection.infon, function(o){return o['@key']=='flatter';})['#text'];
      //
      //       $('#points').html(points);
      //       $('#partner-username').html(partner);
      //       $('#partner-level').html(partner_level);
      //       $('#flatter').html(flatter);
      //       $('#partner-results').slideDown(function() {
      //         update_score();
      //         #<{(| Show underlines |)}>#
      //         _.each(data.collection.document.passage, function(passage, passage_idx) {
      //           YPet[passage_idx].currentView.drawBioC(passage, true);
      //         });
      //       });
      //     },
      //     error: function(error_res) {
      //       #<{(| No response so show them the annotation message |)}>#
      //       _.each(passages, function(passage, passage_idx) {
      //         YPet[passage_idx].currentView.drawBioC(null, true);
      //       });
      //       window.scrollTo(0,0);
      //       $('#new-alert-congrats').slideDown()
      //     }
      //   });
      // };
      //
      // $('#quest-submit').one('click', function(evt) {
      //   #<{(| AJAX submit all of the current Annotaitons |)}>#
      //   var self = this;
      //       ann_dict = {};
      //
      //   #<{(| Iterate over each of the paragraphs or annotatable sections on the page |)}>#
      //   _.each(passages, function(passage, passage_idx) {
      //     ann_dict[+_.find(passage.infon, function(o){return o['@key']=='id';})['#text']] = YPet[passage_idx].currentView.collection.parentDocument.get('annotations').toJSON();
      //   });
      //
      //   #<{(| Add the section to the objects |)}>#
      //   _.each(_.keys(ann_dict), function(section_pk) {
      //     ann_dict[section_pk] = _.map(ann_dict[section_pk], function(obj) { return _.extend(obj, {'section_pk': section_pk}) })
      //   })
      //
      //   #<{(| Do not save data if the annotation is empty |)}>#
      //   annotations = _.flatten(_.values(ann_dict));
      //   annotations = _.difference(annotations, _.where(annotations, {'text': ""}));
      //   annotations = _.map(annotations, function(o) { return _.omit(o, 'opponent');});
      //   annotations =_.map(annotations, function(o) { return _.omit(o, 'words');});
      //
      //   #<{(| Submit Task over ajax, then show correct page (new / gm / partner compare) |)}>#
      //   $.ajax({
      //     type: 'POST',
      //     url: '/task/entity-recognition/quest/{{task.pk}}/{{document.pk}}/submit/',
      //     headers: {'X-CSRFTOKEN': '{{csrf_token}}'},
      //     contentType: "application/json; charset=utf-8",
      //     data:  JSON.stringify(annotations),
      //     dataType: 'json',
      //     cache: false,
      //     async: false,
      //     success: function() {
      //       $('#doc-talk').show();
      //       $(self).hide();
      //       show_results();
      //     }
      //   });
      //
      // });
  }
});

