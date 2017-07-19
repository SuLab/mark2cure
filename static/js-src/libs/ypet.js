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

var channel = Backbone.Radio.channel('ypet');

/*
 *  Models & Collections
 */

NERAnnotationTypeList = Backbone.Collection.extend({
  /* Very simple collection to store the type of
   * Annotations that the application allows
   * for paragraphs */
  model: Backbone.Model.extend({}),
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
});

NERWordList = Backbone.Collection.extend({
  /* Common utils to perform on an array of Word
   * models for house keeping and search */
  model: NERWord,
});

NERAnnotation = Backbone.RelationalModel.extend({
  /* A User or Opponent Annotation (contains NERWordList) */
  defaults: {
    /* An annotation doesn't exist when removed so
     * we can start them all off at 0 and not need to
     * mix in a null type */
    type_id: 0,
    text: '',
    start: null,
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
    var words_len = this.get('words').length;

    if(words_len == 0 && this.get('text') != '') {
      /* If the Annotation was made using text */
      alert('implement text >> words fucntion');

    } else if(words_len > 0 && this.get('text') == '') {
      /* If the Annotation was made using NERWordList */

      var ann = this.sanitizeAnnotation(this.get('words').pluck('text').join(' '), this.get('words').first().get('start'));
      this.set('text', ann.text);
      this.set('start', ann.start);

    } else {
      alert('major issue');
    }

    /* (TODO) Only if user created */
    channel.trigger('ypet:footer:search:set', this.get('text'));

  },

  sanitizeAnnotation: function(full_str, start) {
    /* Return the cleaned string and the (potentially) new start position */
    var str = _.str.clean(full_str).replace(/^[^a-z\d]*|[^a-z\d]*$/gi, '');
    return {'text':str, 'start': start+full_str.indexOf(str)};
  },

  toggleType: function() {
    /* Removes (if only 1 Annotation type) or changes
     * the Annotation type when clicked after existing */
    if( this.get('type_id') == NERAnnotationTypes.length-1 || this.get('text') == '') {
      this.destroy();
    } else {
      this.set('type_id', this.get('type_id')+1 );
    }
  },

  draw: function() {
    var annotation_type = NERAnnotationTypes.at(this.get('type_id'));
    var words_len = this.get('words').length;

    this.get('words').each(function(word, word_index) {
      if(word_index == words_len-1) { word.set('neighbor', true); }
      word.trigger('highlight', annotation_type.get('color'));
    });

    /* Draw all the basic background or underlines */
    // word.trigger('underline', {'color': annotation_type.get('color')});
    // if(annotation.get('opponent')) {
    //   var words = annotation.get('words')
    //   var author_annotations = parent_document.get('annotations');
    //
    //   var anns = []
    //   author_annotations.each(function(main_ann) {
    //     if(main_ann.get('words').contains(words.first()) || main_ann.get('words').contains(words.last())) {
    //       anns.push(main_ann.cid);
    //     }
    //   });
    //
    //   if(_.uniq(anns).length > 1) {
    //     #<{(| 2 Different Parent Annotations |)}>#
    //     words.each(function(word) {
    //       if(word == words.last()) {
    //         word.trigger('underline-space', {'color': '#fff', 'last_word': true});
    //       } else {
    //         word.trigger('underline-space', {'color': annotation_type.get('color'), 'last_word': false});
    //       }
    //     });
    //   }
    // }
  }
});

NERAnnotationList = Backbone.Collection.extend({
  /* Utils for the Paragraph Annotations lists
   * collectively */
  model: NERAnnotation,

  initialize: function(options) {
    /* When a new Annotation has been added */
    this.listenTo(this, 'add', function(annotation) { annotation.draw(); });
    this.listenTo(this, 'change:type_id', function(annotation) { annotation.draw(); });
    this.listenTo(this, 'remove', function(annotation, collection) {
      /* (TODO) Can probably be optimzed rather than a total refresh */
      channel.trigger('ypet:paragraphs:redraw');
    });
  },

  draw: function() {
    this.each(function(annotation) { annotation.draw(); });
  }

});

NERParagraph = Backbone.RelationalModel.extend({
  /* Foundational model for tracking everything going
   * on within a Paragraph like Words and Annotations */
  defaults: {
    text: '',
    offset: 0,
  },

  relations: [{
    type: 'HasMany',
    key: 'words',
    relatedModel: NERWord,
    collectionType: NERWordList,
  }, {
    type: 'HasMany',
    key: 'annotations',
    relatedModel: NERAnnotation,
    collectionType: NERAnnotationList,
    collectionOptions: {foo: 6969},
  }, {
    type: 'HasMany',
    key: 'opponent_annotations',
    relatedModel: NERAnnotation,
    collectionType: NERAnnotationList,
    collectionOptions: {foo: 123},
  }],

  initialize : function() {
    // console.log('ner para init', this.get('opponent_annotations') );
    //
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

    this.set('offset', +this.get('offset'));

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
});

NERDocumentResult = Backbone.RelationalModel.extend({
  /* Provides the information for comparison between
   * the User and a selected opponent */
  defaults: {
    'flatter': "",
    'award': {
      'amount': 0
    },
    'opponent': null
  }
})

NERQuestResult = Backbone.RelationalModel.extend({
  /* Provides the information for Quest results */
  defaults: {
    'task': {
      'pk': null,
      'name': ''
    },
    'group': {
        'pk': null,
        'name': '',
        'stub': '',
        'enabled': false,
        'description': ''
    },
    'award': {
      'pk': null,
      'amount': 0
    },
    'uqr_created': false,
    'uqr_pk': null
  }
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
    active: false, // If is the currently addressed Document in the Quest
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

  url: function() {
    return '/api/ner/quest/'+ this.quest_pk +'/';
  },

  initialize: function(options) {
    this.quest_pk = options.quest_pk;


    this.listenTo(this, 'add', function(model) {
      // (TODO) Cleanup any attributes?
    });

    // this.listenTo(this, 'sync', function() {
    //   if( this.is_quest_complete() ) {
    //     // Jump to review page b/c the Quest is complete!
    //     this.trigger('quest:complete');
    //   }
    // });

  },

  is_quest_complete: function() {
    var quest_completed = this.pluck('quest_completed').every(function(v){ return v == 1; }),
        all_documents_completed = this.pluck('document_view_completed').every(function(v){ return v == 1; });
    return quest_completed == false && all_documents_completed;
  },

  get_active: function() {
    var available = this.where({'document_view_completed': 0});
    if(available.length == 0) {

      if(this.is_quest_complete()) {
        this.trigger('quest:complete');
      } else {
        this.trigger('quest:finished');
      }

      return null;
    }

    var m = available[_.random(0, available.length-1)]
    m.set('active', true);
    return m;
  }
})

/*
 * Views
 */

NERWordView = Backbone.Marionette.View.extend({
  /* View for all direct actions on a word
  * - Model = NERWord
  * - Collection = None
  * - options.section_pk = PK of current section (passage / paragraph)
  */
  template: _.template('<% if(neighbor) { %><%= text %><% } else { %><%= text %> <% } %>'),
  tagName: 'span',

  modelEvents: {
    'change:neighbor': function() { this.render(); },
    'change:disabled': 'onChangeDisabled',
    'change:latest': 'onChangeLatest',

    'unclick': 'onUnClick',
    'highlight': 'onHighlight',
    'underline': 'onUnderline',
    'underline-space': 'onUnderlineSpace'
  },

  /* These events are only triggered when over
   * a span in the paragraph */
  events : {
    'mousedown' : 'onMouseDown',
    'mouseover' : 'onMouseOver',
    'mouseup'   : 'onMouseUp',
  },

  /* Triggers the proper class assignment
   * when the word <span> is redrawn */
  onRender: function() {
    this.$el.css(this.model.get('neighbor') ?
      {'margin-right': '5px', 'padding-right': '0'} :
      {'margin-right': '0px', 'padding-right': '.25rem'});
  },

  onChangeDisabled: function() {
    /* Apply not-allowed CSS to the $el */
    if(this.model.get('disabled')) {
      this.$el.css('cursor', 'not-allowed');
    }
  },

  onChangeLatest: function(model, value, options) {
    /* When the latest attribute is set, means they're hover over so highlight temporarily */
    if(this.model.get('latest')) {
      this.model.trigger('highlight', '#D1F3FF');
    }
    if(options.force) {
      this.model.trigger('highlight', '#fff');
    }
  },

  onUnClick: function() {
    /* Set $el's background-color to white, then trigger mouseup
     * Ex. for when we want to force remove a selection (during training) */
    var $el = this.$el;
    $el.animate({backgroundColor: '#fff'}, 750, function() {
      $el.trigger('mousedown').trigger('mouseup');
    });
  },

  onHighlight: function(color) {
    /* Apply background-color CSS to the $el */
    this.$el.css({'backgroundColor': color});
  },

  onUnderline: function() {
    /* (UNCLEAR) create a new absolute positioned colored div */
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

  onUnderlineSpace: function() {
    /* (UNCLEAR) create a new absolute positioned white div */
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

  onMouseDown: function(evt) {
    /* Set the Word's latest attribute to 1 */
    if(this.model.get('disabled')) { return; };
    this.model.set({'latest': 1});
  },

  onMouseOver: function(evt) {
    /* Set array of Word's latest attribute to Date.now() to track drag */

    if(this.model.get('disabled')) { return; };
    var word_model = this.model,
        word_collection = this.model.collection;

    /* You're dragging if another word has a latest timestamp */
    if(_.compact(word_collection.pluck('latest')).length) {
      if(_.isNull(word_model.get('latest'))) { word_model.set({'latest': Date.now()}); }

      /* If the hover doesn't proceed in ordered fashion
       * we need to "fill in the blanks" between the words */

      // Index position of the currently hovered word
      var current_word_idx = word_collection.indexOf(word_model);
      // Index position of the first word that was selected
      var first_word_idx = word_collection.indexOf( word_collection.find(function(w) { return w.get('latest') == 1; }) );

      /* Select everything from the starting to the end without
       * updating the timestamp on the first_word */
      var starting_positions = first_word_idx <= current_word_idx ? [first_word_idx, current_word_idx+1] : [first_word_idx+1, current_word_idx];
      var selection_indexes = _.range(_.min(starting_positions), _.max(starting_positions));
      _.each(_.without(selection_indexes, first_word_idx), function(idx) { word_collection.at(idx).set('latest', Date.now()); });

      /* If there are extra word selections up or downstream
       * from the current selection, remove those */
      var last_selection_indexes = _.map(word_collection.reject(function(word) { return _.isNull(word.get('latest')); }), function(word) { return word_collection.indexOf(word); });
      var remove_indexes = _.difference(last_selection_indexes, selection_indexes);

      var remove_word, ann;
      _.each(remove_indexes, function(idx) {
        remove_word = word_collection.at(idx);
        remove_word.set('latest', null, {'force': true});

        // (TODO) What is going on here?
        // ann = word.get('parentAnnotation');
        // if(ann) { ann.collection.drawAnnotation(ann); }

      });
    }
  },

  onMouseUp: function(evt) {
    /* After click, drag, etc send words to AnnotationList as Annotation */

    if(this.model.get('disabled')) { return; };

    var selected = this.model.collection.filter(function(w) { return w.get('latest') });
    channel.trigger('ypet:paragraph:set:annotations', {
      'section_pk': this.getOption('section_pk'),
      'words': selected,
      'word': this.model
    });

    // Reset all latest attributes for future selections
    this.model.collection.each(function(w) { w.set('latest', null); });
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
  childViewEventPrefix: 'word',
  childViewOptions: function() {
    return {'section_pk': this.getOption('section_pk')}
  }
});

NERParagraphView = Backbone.Marionette.View.extend({
  /* Container item for the individual words
   * this.model = NERParagraph
   * this.collection = None
  */
  template: _.template('<div class="paragraph-content p-3"></div>'),
  className: function() {
    var enabled_status = this.getOption('enabled') ? 'ypet-enabled' : 'ypet-disabled';
    return ['paragraph-box', 'my-3', enabled_status].join(' ');
  },

  regions: {
    'words': '.paragraph-content'
  },

  initialize: function() {
    this.listenTo(channel, 'ypet:paragraphs:redraw', this.render);
  },

  onRender : function() {
    if(!this.options.enabled) {
      /* If you're showing a partner's results, disallow highlighting */
      // this.$el.css({'color': '#000', 'cursor': 'default'});;
      this.model.get('words').each(function(w) {
        w.set('disabled', true);
      });
      this.$el.css('cursor', 'not-allowed');
    }

    this.showChildView('words', new NERWordsView(
      { 'collection': this.model.get('words'),
        'enabled': this.options.enabled,
        'section_pk': this.model.get('pk')
      }
    ));

    this.model.get('annotations').draw();
    this.model.get('opponent_annotations').draw();

    // var has_opponent =  _.contains(this.model.get('annotations').pluck('self'), false);
    // if (has_opponent) {
    //   #<{(| Warning: this allows paragraph specific disabiling |)}>#
    //   this.triggerMethod('view:only');
    // }
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
    if(this.getOption('mode') == 'er') { enabled = false; }
    if(this.getOption('mode') == 'ner' && this.getOption('review') == true) { enabled = false; }
    return {
      'enabled': enabled
    }
  },
  initialize: function() {
    if(this.model) {
      this.collection = this.model.get('passages');
    } else {
      this.collection = new NERParagraphList({});
    }

    channel.reply('ypet:paragraph:user:annotations', this.getUserAnnotations, this);
    channel.reply('ypet:paragraph:opponent:annotations', this.getOpponentAnnotations, this);

    this.listenTo(channel, 'ypet:paragraph:set:annotations', function(obj) {
      var passage = this.model.get('passages').findWhere({'pk': obj.section_pk});

      if(obj.words.length == 1) {

        passage.get('annotations').each(function(annotation) {
          if(annotation.get('words').contains(obj.word)) {
            annotation.toggleType();
          };
        });

      } else {
        passage.get('annotations').create({'words': obj.words});
      }
      // passage.get('annotations').each(function(annotation) {
      //   console.log( annotation.get('words').contains(obj.word) );
      // });

      // if(selected.length == 1 && word_model.get('parentAnnotation') ) {
      //   word_model.get('parentAnnotation').toggleType();
      // } else {
      //   #<{(| if selection includes an annotation, delete that one |)}>#
      //   _.each(selected, function(w) {
      //     if(w.get('parentAnnotation')) {
      //       w.get('parentAnnotation').destroy();
      //     }
      //   });
      //   #<{(| Take the selected words as a whole and make a new Annotation with them |)}>#
      //   word_model.get('parentDocument').get('annotations').create({'words': selected});
      // };

    });

  },

  getUserAnnotations: function(section_pk) {
    /* Return the Paragraph's Annotations */
    var passage = this.model.get('passages').findWhere({'pk': section_pk});
    return passage.get('annotations')
  },

  getOpponentAnnotations: function(section_pk) {
    /* Return the Paragraph's Opponent Annotations */
    var passage = this.model.get('passages').findWhere({'pk': section_pk});
    return passage.get('opponent_annotations')
  },
});

NERDocumentResultsView = Backbone.Marionette.View.extend({
  /* ...
   * this.model = NERDocumentResult
   * this.collection = None
   */
  template: '#ypet-document-results-template',
  className: 'row',

  initialize: function() {
    var document_pk = this.model.get('pk');
    this.model = new NERDocumentResult({});

    /* Replace the model. NERDocument >> NERDocumentResult */
    var self = this;

    $.ajax({
      url:'/task/entity-recognition/'+ this.options.task_pk +'/'+ document_pk +'/results.json',
      dataType: 'json',
      headers: {'X-CSRFTOKEN': this.options.csrf_token},
      success: function(data) {
        self.model = new NERDocumentResult(data);
        channel.trigger('ypet:navigation:score:update');
        self.render();
      },
      error: function(error_res) {
        /* (TODO) Show Error */
        window.scrollTo(0,0);
      }
    });
  },

  onRender: function() {
    window.scrollTo(0,0);
  }
});

NERQuestCompletedView = Backbone.Marionette.View.extend({
  /* ...
   * this.model = NERQuestResult
   * this.collection = None
   */
  template: '#ypet-quest-completed-template',
  className: 'row',

  ui: {
    'next': '#next-quest'
  },

  events: {
    'mousedown @ui.next': 'nextQuest'
  },

  initialize: function() {
    this.model = new NERQuestResult({});
    var self = this;

    // This is either creating for the first time or just posting to fetch the UQR info
    $.ajax({
      type: 'POST',
      url: '/task/entity-recognition/quest/'+this.options.task_pk+'/submit/',
      headers: {'X-CSRFTOKEN': this.options.csrf_token},
      success: function(data) {
        self.model = new NERQuestResult(data)
        channel.trigger('ypet:navigation:score:update');
        self.render();
      }
    });
  },

  nextQuest: function() {
    /* Use API to determine next Quest for User */

    $.ajax({
      type: 'GET',
      url: '/api/ner/list/'+ this.model.get('group').pk +'/quests/',
      headers: {'X-CSRFTOKEN': this.options.csrf_token},
      success: function(data) {
        var set = _.filter(data, function(q) {
          return q.progress.completed == false && q.user.completed == false && q.user.enabled == true;
        });

        if(set.length) {
          window.location = '/task/entity-recognition/quest/' + set[0].id;
        } else { window.location = '/dashboard/'; }

      }
    });

  }
});

NERLoadingView = Backbone.Marionette.View.extend({
  /* Initial HTML before a REExtractionList is available
  * - Model = None
  * - Collection = None
  */
  template: '<p>Text Loading...</p>'
});

NERProgressItem = Backbone.Marionette.View.extend({
  /* Drop down list of REChoices */
  template: _.template('&#8226;'),
  tagName: 'li',
  className: 'list-inline-item',

  onRender: function() {
    if ( this.model.get('document_view_completed') ) {
      this.$el.addClass('completed');
    }
    if ( this.model.get('active') ) {
      this.$el.addClass('active');
    }
  }
});

NERProgressView = Backbone.Marionette.CollectionView.extend({
  /* Parent list for REProgressItem */
  tagName: 'ul',
  className: 'list-unstyled list-inline',
  childView: NERProgressItem,
  childViewEventPrefix: 'progress',
});

NERNavigationView = Backbone.Marionette.View.extend({
  /* The Progress indicator view for all interactions on the Quest
   * this.model = None
   * this.collection = NERDocumentList (Quest Documents)
   */
  template: '#ypet-navigation-template',
  className: 'row justify-content-between',

  regions: {
    'progress': '#progress-bar',
  },

  ui: {
    'score': '#score'
  },

  initialize: function() {
    this.listenTo(channel, 'ypet:navigation:score:update', this.updateScore);
  },

  onRender: function() {
    var self = this;

    // this.od = new Odometer({
    //   el: self.ui.score.first(),
    //   value: 0,
    //   format: '(,ddd)',
    //   theme: 'minimal'
    // });

    this.showChildView('progress', new NERProgressView({'collection': this.collection}));
    this.updateScore();
  },

  updateScore: function() {
    var self = this;

    $.ajax({
      type: 'GET',
      url: '/u/points/',
      headers: {'X-CSRFTOKEN': this.options.csrf_token},
      cache: false,
      async: false,
      success: function(data) {
        self.ui.score.html(data.points);
        // self.od.update(data.points);
      }
    });

  }
});

NERFooterHelpView = Backbone.Marionette.View.extend({
  /* The list of links for getting task instructions
   * this.model = None
   * this.collection = None
   */
  template: '#ypet-footer-help-template',
  templateContext: function() {
    return {'mode': this.options.mode};
  },
  className: 'row'
})

NERFooterConfirmView = Backbone.Marionette.View.extend({
  /* The Actionable button that submits the document
   * this.model = None
   * this.collection = None
   */
  template: _.template('Submit'),
  className: 'btn btn-block confirm-button',
  tagName: 'p',

  triggers: {
    'mousedown': 'confirm:submit'
  }
});

NERFooterNextView = Backbone.Marionette.View.extend({
  /* The Actionable button that reviews the document discussion
   * page or goes to the next avilable NER Document
   * this.model = None
   * this.collection = None
   */
  template: '#ypet-footer-next-template',
  className: 'row',

  ui: {
    'next': '.next-doc'
  },

  triggers: {
    'mousedown @ui.next': 'next'
  }
});

NERFooterSearchView = Backbone.Marionette.View.extend({
  /* The link that allows people to look up a term independently
   * this.model = None
   * this.collection = None
   */
  template: '#ypet-footer-search-template',
  className: 'row',

  ui: {
    'link': 'a',
    'small': 'small'
  },

  initialize: function() {
    this.listenTo(channel, 'ypet:footer:search:set', this.showText);
  },

  onRender: function() {
    this.$el.hide();
  },

  showText: function(text) {
    this.$el.show();
    var url = 'https://www.google.com/search?q='+text;
    this.ui.link.attr('href', url);
    this.ui.small.text(_.str.truncate(text, 36));
  }
});

NERFooterView = Backbone.Marionette.View.extend({
  template: '#ypet-footer-template',
  className: 'row my-3 justify-content-between',
  childViewEventPrefix: 'ner:footer',

  regions: {
    'help': '#ypet-footer-help',
    'confirm': '#ypet-footer-confirm',
    'search': '#ypet-footer-search'
  },

  initialize: function() {
    /* Is this a submit or next state
    * States:
    * review - confirm (still viewing a doc)
    * review - next (done reviewing a doc)
    * */
  },

  onRender: function() {
    this.showChildView('help', new NERFooterHelpView({'mode': this.options.mode}));


    this.showChildView('confirm', new NERFooterConfirmView({}));
    if(this.getOption('review') == true){
      this.showChildView('confirm', new NERFooterNextView({}));
    }

    // (TODO) Pass NERAnnotation item into this
    this.showChildView('search', new NERFooterSearchView({'model': false}));
  }
});

YPet = Backbone.Marionette.View.extend({
  /* The top level view for all interactions on text
   * this.model = active NERDocument being annotated
   * this.collection = NERDocumentList (Quest Documents)
   */
  template: '#ypet-template',
  className: 'row justify-content-center',

  regions: {
    'navigation': '#ypet-navigation',
    'results': '#ypet-results',
    'text': '#ypet-text',
    'footer': '#ypet-footer'
  },

  initialize: function() {
    if(!this.options.training && !this.collection) {
      this.collection = new NERDocumentList({'quest_pk': this.options.task_pk});
      this.collection.fetch();
    }

    this.listenTo(this.collection, 'sync', function() {
      this.model = this.collection.get_active();
      if(this.model) { this.render(); }
    });
  },

  collectionEvents: {
    'quest:complete': function() {
      // When the user just completed the last NER Task
      this.getRegion('text').empty();
      this.getRegion('footer').empty();
      this.showChildView('results', new NERQuestCompletedView(this.options));
    },
    'quest:finished': function() {
      // When the user revists a previously completed Quest
      this.getRegion('text').empty();
      this.getRegion('footer').empty();
      this.showChildView('results', new NERQuestCompletedView(this.options));
    }
  },

  childViewEvents: {
    'ner:footer:next': function() {
      /* Done reviewing the current document, go to the next */
      this.model = this.collection.get_active();
      this.options['review'] = false
      this.render();
    },
    'ner:footer:confirm:submit': function() {
      /* Submit the User's annotations for the current Document */
      var self = this;
          ann_dict = {};

      /* Iterate over each of the paragraphs or annotatable sections on the page */
      _.each(this.model.get('passages').pluck('pk'), function(section_pk) {
        ann_dict[section_pk] = channel.request('ypet:paragraph:user:annotations', section_pk).toJSON();
        ann_dict[section_pk] = _.map(ann_dict[section_pk], function(obj) { return _.extend(obj, {'section_pk': +section_pk}) })
      });

      /* Do not save data if the annotation is empty */
      annotations = _.flatten(_.values(ann_dict));
      annotations = _.difference(annotations, _.where(annotations, {'text': ""}));
      annotations = _.map(annotations, function(o) { return _.omit(o, 'opponent');});
      annotations =_.map(annotations, function(o) { return _.omit(o, 'words');});

      /* Submit Task over ajax, then show correct page (new / gm / partner compare) */
      $.ajax({
        type: 'POST',
        url: '/task/entity-recognition/quest/'+this.options.task_pk+'/'+this.model.get('pk')+'/submit/',
        headers: {'X-CSRFTOKEN': this.options.csrf_token},
        contentType: "application/json; charset=utf-8",
        data:  JSON.stringify(annotations),
        dataType: 'json',
        cache: false,
        async: false,
        success: function() {
          self.options['review'] = true
          self.render();
        }
      });

    }
  },

  onRender: function() {
    if(this.collection) {
      this.options['model'] = this.model;
      this.options['collection'] = this.collection;
      this.showChildView('text', new NERParagraphsView(this.options));

      if(this.options.mode == 'ner') {
        this.showChildView('navigation', new NERNavigationView(this.options));
        this.showChildView('footer', new NERFooterView(this.options));

        if(this.getOption('review') == true) {
          this.showChildView('results', new NERDocumentResultsView(this.options));
        }

      } else if (this.options.mode == 're') {
        // var concept_uids = [this.options.concepts['c1'].id, this.options.concepts['c2'].id];
      }
    } else {
      this.showChildView('text', new NERLoadingView());
    };
  }
});

